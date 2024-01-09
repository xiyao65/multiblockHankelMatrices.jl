# this file is to make the multi_level Hankel matrix explicitly
# module multihankel
# import Pkg;
# Pkg.add("DSP")
# Pkg.add("FFTW")
# Pkg.add("AbstractFFTs")
using DSP
using BenchmarkTools
using FFTW
using FFTW: Plan
using Random
# include("multiblockHankel.jl")

using Random

import Base: SubArray, Colon, to_indices,hcat,vcat,muladd,prod!,convert,size,getindex,adjoint,conj
import Base: adjoint, convert, transpose, size, similar, copy, getproperty, inv, sqrt, copyto!, reverse, conj, zero, fill!, checkbounds, real, imag, isfinite, DimsInteger, iszero
import LinearAlgebra: Cholesky, DimensionMismatch, cholesky, cholesky!, eigvals, inv, ldiv!,
    mul!, pinv, rmul!, tril, triu,norm,normalize,svd,Diagonal
using LinearAlgebra: LinearAlgebra, Adjoint, Factorization, factorize,qr
import Base: ==, +, -, *, \

export Hankel

# # # # # # # # # # # # # # # # # # #   for FFT  # # # # # # # # # # # #
function irlblr(H::AbstractHankel{T}, L::Int64,tol::Float64,maxit::Int64,kn::Int64,U0::AbstractMatrix,V0::AbstractMatrix) where T

   
   ##prepare for FFTW
   

   
   
    n_d=H.n_d
    p_d=H.p_d
    q_d=n_d.-p_d.+1
      # hankelmatrix is pp*qq dimension
    nu     = L # the size of extreme singular value
    m      = prod(p_d)
    n      = prod(q_d)
    # m_b    = min(nu+5, 2*nu, n)  # Working dimension size
    m_b    = kn # Working dimension size
    mprod  = 0
    it     = 0  # restart lanczos iterations
    j      = 1  
    k      = nu
    smax   = 1 # check for convergence
    fn =1
   
    V  = zeros(ComplexF64,n,m_b)
    U  = zeros(ComplexF64,m,m_b)
    v  = zeros(T,n) # F/p->v s->u
    B  = zeros(Float64,m_b,m_b)
    error_end=10
    V[:,1]  =normalize(rand(T,n)) # Initial vector
   @inbounds begin
    while(it < maxit)
      if(it>0) j=k end  
      # U[:,j] = matrixvector.stringvector(n_d,S,V[:,j],dft_s,idft_s)+U0*(V0'*V[:,j])
      U[:,j] = H*V[:,j]+U0*(V0'*V[:,j])
      mprod+=1
      mprod+=1
      mprod+=1
      if(it>0)
        U[:,j] = orthog(U[:,1:j-1],U[:,j])
      end
       u = norm(U[:,j])
      uinv = invcheck(u)
      U[:,j] = uinv*U[:,j]
      # Lanczos process
          while(j<m_b)
           # v = matrixvector.adstringvector(n_d,S,U[:,j],dft_s,idft_s)+V0*(U0'*U[:,j])
           v = H'*U[:,j]+V0*(U0'*U[:,j])
            mprod+=1
            v = v - u*V[:,j]
            v = orthog(V[:,1:j],v)
            fn = norm(v)
            fninv= invcheck(fn)
            v = fninv * v
                if(j<m_b-1)
                  V[:,j+1] = v
                  B[j,j] = u
                  B[j,j+1] = fn
                #  U[:,j+1] =  matrixvector.stringvector(n_d,S,V[:,j+1],dft_s,idft_s)+U0*(V0'*V[:,j+1])
                 U[:,j+1] = H*V[:,j+1]+U0*(V0'*V[:,j+1])
                  mprod+=1
                  U[:,j+1] = U[:,j+1] - fn*U[:,j]
                  # ...with full reorthogonalization
                  U[:,j+1] = orthog(U[:,1:j],U[:,j+1])
                  u = norm(U[:,j+1])
                  uinv = invcheck(u)
                  U[:,j+1] = uinv * U[:,j+1]
                else
                  B[j,j] = u
   
               end
          j+=1
          end
          # end of one iteration of lanczos
         F = svd(B)
         R = fn*F.U[m_b-1,:] # Residuals
         error_end=R[nu-1]
         if(it<1)
           smax = F.S[1]
         else
           smax = max(F.S[1],smax)
         end
         # conv = count( <(tol*smax),abs.(R[1:nu]))
         conv = count( <(tol),abs.(R[1:nu]))
         # conv = k'
         if (conv<nu)
           # k = max(conv+1, k)
           k = max(conv+nu, k)
           k = min(k, m_b-3)
           k =max(k,1)
         else
           break
         end
         # resart with this part
         V[:,1:k-1] = V[:,1:m_b]*F.V[:,1:k-1]
         V[:,k] = v
         B[1:k-1,1:k-1]=Diagonal(F.S[1:k-1])
         B[1:k,k]=R[1:k]
         U[:,1:k]=U[:,1:m_b]*F.U[:,1:k]
         it+=1
       end #(while end it )
   if (it<1)
    F = svd(B)#(first svd B)
   end
    U = U[:,1:m_b]*F.U[:,1:nu]
    V = V[:,1:m_b]*F.V[:,1:nu]
   
    # println(["#mprod"              string(mprod)
    #            "#it"              string(it)
    #            "#error"              string(error_end)
    #            ])
    return (U,F.S[1:nu],V,it, mprod)
   
   end # for inbounds
   
   end

   function pad(x::AbstractArray{T},n::Vector{Int},x_n::Vector{Int}) where T
    nss=zeros(T,tuple(n...))
    x_array=reshape(reverse(x),tuple(x_n...))
    in=axes(x_array)
    nss[in...] =x_array
    result=vec(nss)
    return result
end
function riemannsvd(U::Array{<:Number},V::Array{<:Number},H::AbstractHankel{T}) where T
    r=size(U,2)
    HV=H*V
    adHU=H'*U
    G=adjoint(U)*HV
    m1=adjoint(V)*adHU
    B=adHU-V*m1
    m2=adjoint(U)*HV
    C=HV-U*m2
    q_1,r_1=qr(B)
    q_2,r_2=qr(C)
    Middle=ComplexF64.(zeros(T,2*r,2*r))
    Middle[1:r,1:r]=G
    Middle[(r+1):end,1:r]=r_2
    Middle[1:r,(r+1):end]=adjoint(r_1)
    U_M,sigma_M,V_M=svd(Middle)
    U=(U*U_M[1:r,:]+q_2*U_M[r+1:end,:])[:,1:r]
    V=(V*V_M[1:r,:]+q_1*V_M[r+1:end,:])[:,1:r]
    S=sigma_M[1:r]
    return (U,S,V)
end
function extract(x::AbstractArray{T},n::Vector{Int},x_n::Vector{Int}) where T
    x_array=reshape(reverse(x),tuple(n...))
    tmp=zeros(T,tuple(x_n...))
    in=axes(tmp)
    tmp=x_array[in...]
    result=reverse(vec(tmp))
    return result
end
function stringvector(s::Vector{T},p_s::Vector{Int64},q_s::Vector{Int64},x::AbstractArray{<:Number},y::AbstractArray{<:Number},α::Number) where T <:Number
  # function stringvector(s::Vector{T},p_s::Vector{Int64},q_s::Vector{Int64},x::Vector{T},y::Vector{T},α::Number) where T
    # stringvector(A.element,A.p_d,A.q_d,x,y,α )
    #    s=ComplexF64.(s)
       @inbounds for j in 1:prod(q_s)
           tmp = α * x[j]
       for i in 1:prod(p_s)
           y[i] =muladd(tmp,s[i+j-1],y[i])
       end
       end
       return y
end
function fstringvector(tmp::Array{T},n_d::Vector{Int64},p_d::Vector{Int64},q_d::Vector{Int64},
    x::AbstractArray,y::AbstractArray,α::Number,dft_s::Plan{ComplexF64},idft_s::Plan{ComplexF64}) where T

  # fstringvector(B.tmp,B.n_d,B.p_d,B.q_d,x,y,α,B.dff_s,B.idft_s)


    tmp1=pad(x,n_d,q_d)
    y_hat=idft_s*(dft_s*tmp.*(dft_s*tmp1))
    y_hat=extract(y_hat,n_d,p_d)
    y= muladd(α, y_hat, y)
    y=maybereal(T,y)
    return y
end

   # # # # # # # # # # # # # # # # # # #   for FFT  # # # # # # # # # # # #

abstract type AbstractHankel{T<:Number} <: AbstractMatrix{T} end


size(A::AbstractHankel) = (size(A, 1), size(A, 2))


function getindex(A::AbstractHankel, i::Integer)
    return A[mod(i - 1, size(A, 1)) + 1, div(i - 1, size(A, 1)) + 1]
end

convert(::Type{AbstractMatrix{T}}, S::AbstractHankel) where {T} = convert(AbstractHankel{T}, S)
convert(::Type{AbstractArray{T}}, S::AbstractHankel) where {T} = convert(AbstractHankel{T}, S)




# Fast application of a general Hankel matrix to a column vector via FFT
function mul!(y::StridedVector, A::AbstractHankel, x::StridedVector, α::Number, β::Number)
    m, n = size(A)
    if length(y) != m
        throw(DimensionMismatch(
            "first dimension of A, $(m), does not match length of y, $(length(y))"
        ))
    end
    if length(x) != n
        throw(DimensionMismatch(
            "second dimension of A, $(n), does not match length of x, $(length(x))"
        ))
    end
    if iszero(β)
        fill!(y, 0)
    else
        rmul!(y, β)
    end
   
    N = m + n - 1
    if N < 512   && size(A.p_d,1)==1
         # Small case: don't use FFT
        y=stringvector(A.element,A.p_d,A.q_d,x,y,α )

    else
        # Large case: use FFT
        y=fstringvector(vec(A.element),A.n_d,A.p_d,A.q_d,x,y,α,A.dft_s,A.idft_s)
    end
    return y
end


function mul!(C::StridedMatrix, A::AbstractHankel, B::StridedMatrix, α::Number, β::Number)
    l = size(B, 2)
    if size(C, 2) != l
        throw(DimensionMismatch("input and output matrices must have same number of columns"))
    end
    for j = 1:l
    C[:,j]   = mul!(view(C, :, j), A, view(B, :, j), α, β)
    end
    return C
end

# General Hankel matrix
"""
    Hankel

A Hankel matrix.
"""
struct Hankel{T<:Number} <: AbstractHankel{T}
    element::Array{T}
    n_d::Vector{Int64}
    p_d::Vector{Int64}
    q_d::Vector{Int64}
    dft_s::Plan
    idft_s::Plan

    function Hankel{T}(element::Array{T}, p_d::Vector{Int64}) where {T<:Number}
        # if first(vc) != first(vr)
        #     error("First element of the vectors must be the same")
        # end
        n_d=collect(size(element))
        q_d=n_d.-p_d.+1
        tmp = vec(element)
        dft_s=plan_fft(tmp)
        idft_s=plan_ifft(tmp)
        return new{T}(element,n_d,p_d,q_d,dft_s,idft_s)
    end
end




"""
    Hankel(vc::AbstractVector, vr::AbstractVector)

Create a `Hankel` matrix from its element and p_d dimension.
"""
function Hankel(element::AbstractArray, p_d::AbstractVector)
    return Hankel{Base.eltype(element)}(element, p_d)
end
function Hankel{T}(element::AbstractVector, p_d::AbstractVector) where {T<:Number}
    return Hankel{T}(convert(Vector{T}, element), convert(Vector{Int64}, p_d))
end


convert(::Type{AbstractHankel{T}}, A::Hankel) where {T} = convert(Hankel{T}, A)
convert(::Type{Hankel{T}}, A::Hankel) where {T} = Hankel(convert(Array{T}, A.element),
                                                               convert(Vector{Int64}, A.p_d))

adjoint(A::Hankel) = Hankel(conj(A.element), A.q_d)
adjoint(A::Hankel{<:Real}) = Hankel(A.element, A.q_d)
transpose(A::Hankel) = Hankel(A.element, A.q_d)
# # Retrieve an entry
function getindex(A::Hankel, i::Integer )
   return vec(A.element)[i]
end

# Size of a general Hankel matrix
function size(A::Hankel, dim::Int)
    if dim == 1
        return prod(A.p_d)
    elseif dim == 2
        return prod(A.q_d)
    elseif dim > 2
        return 1
    else
        error("arraysize: dimension out of range")
    end
end


for fun in (:zero, :conj, :copy, :-, :real, :imag)
    # @eval $fun(A::AbstractHankel)=Hankel($fun(A.element),$fun(A.n_d),$fun(A.p_d),$fun(A.q_d),$fun(A.dft_s),$fun(A.idft_s))
    @eval $fun(A::AbstractHankel)=Hankel($fun(A.element),A.p_d)
end
for op in (:+, :-)
    @eval $op(A::AbstractHankel,B::AbstractHankel)=Hankel($op(A.element,B.element),A.p_d)
end
# # # # # # # # # # # # # # # # # # #   to construct a Hankel matrix explicitly using recursive and it is expensive # # # # # # # # # # # #
function mulHankel(s:: AbstractArray{T},p:: Vector{Int64}) where T
      n_h=size(s)
      k_h=size(n_h,1) # k-dimensional n_1,n_2,n_k array
      p_h=Tuple(p[1:k_h])# q=n_h.-q.+1
      q_h=n_h.-p_h.+1

      if k_h==1
        result = zeros(T,p_h[k_h],q_h[k_h])
        for i in 1:p_h[k_h], j in 1: q_h[k_h]
            result[i,j] = s[i+j-1]
        end
        return result
      else
        result=Matrix{T}
        row_first=Vector{Matrix{T}}(undef,p_h[k_h])
        for i in 1:p_h[k_h]
          row_first[i]=mulHankel(selectdim(s,k_h,i),p[1:k_h-1])
          for j in 1:q_h[k_h]-1
           row_first[i]=hcat(row_first[i],mulHankel(selectdim(s,k_h,i+j),p[1:k_h-1]))
          end
          result=vcat(result,row_first[i])
        end
        result=result[2:end,1:end]
        k_h=k_h-1
    end
    result=convert(Matrix{T},result)
    return result

end

function fullHankel(A:: AbstractHankel{T}) where T
    result=mulHankel(A.element,A.p_d)
    return result
end
# # # # # # # # # # # # # # # # # # #   to construct a Hankel matrix explicitly using recursive and it is expensive # # # # # # #
"""
    maybereal(::Type{T}, x)

Return real-valued part of `x` if `T` is a type of a real number, and `x` otherwise.
"""
maybereal(::Type, x) = x
maybereal(::Type{<:Real}, x) = real(x)

# # # # # # # # # # # # # # # # # # #   truncated svd for multilevel Hankel # # # # # # #
function invcheck(x::Float64)
    eps2  = 2*Base.eps()
    if(x>eps2)
      x = 1/x
    else
      x = 0
     print("Ill-conditioning encountered, result accuracy may be poor")
    end
    return(x)
end

function orthog(X::Matrix{T},Y::Array{T}) where T
      h=zeros(T,size(X,2))
      mul!(h, adjoint(X), Y)
      mul!(Y, X, h, -one(T), one(T))
      return Y
end
# function irlb(S::Array{T}, L::Int64,tol::Float64,maxit::Int64,kn::Int64) where T
function irlb(H::AbstractHankel{T}, L::Int64,tol::Float64,maxit::Int64,kn::Int64) where T
    #  """Estimate a few of the largest singular values and corresponding singular
    # vectors of matrix using the implicitly restarted Lanczos bidiagonalization
    # method of Baglama and Reichel, see:
    # Augmented Implicitly Restarted Lanczos Bidiagonalization Methods,
    # J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005
    # Keyword arguments:
    # tol   -- An estimation tolerance. Smaller means more accurate estimates.
    # maxit -- Maximum number of Lanczos iterations allowed.
    # Given an input matrix A of dimension j * k, and an input desired number
    # of singular values L, the function returns a tuple X with five entries:
    # X[0] A j * nu matrix of estimated left singular vectors.
    # X[1] A vector of length nu of estimated singular values.
    # X[2] A k * nu matrix of estimated right singular vectors.
    # X[3] The number of Lanczos iterations run.
    # X[4] The number of matrix-vector products run.
    # The algorithm estimates the truncated singular value decomposition:
    # A.dot(X[2]) = X[0]*X[1].
    # """
   
   ##prepare for FFTW
   
    S=H.element
    dft_s=H.dft_s
    idft_s=H.idft_s
    # tmp=vec(S)
    # dft_s=plan_fft(tmp)
    # idft_s=plan_ifft(tmp)
   
   
    n_d=H.n_d
    p_d=H.p_d
    q_d=H.q_d
    # n_d=collect(size(S))
    # p_d=fld.(n_d,2).+1
    # q_d=n_d.-p_d.+1
      # hankelmatrix is pp*qq dimension
    nu     = L # the size of extreme singular value
    m      = prod(p_d)
    n      = prod(q_d)
    # m_b    = min(nu+5, 2*nu, n)  # Working dimension size
    m_b    = kn # Working dimension size
    mprod  = 0
    it     = 0  # restart lanczos iterations
    j      = 1  #???? be careful of index
    k      = nu
    smax   = 1 # check for convergence
    fn =1
   
    V  = zeros(ComplexF64,n,m_b)
    U  = zeros(ComplexF64,m,m_b)
    v  = zeros(T,n) # F/p->v s->u
    B  = zeros(Float64,m_b,m_b)
    error_end=10
    V[:,1]  =normalize(rand(T,n)) # Initial vector
   @inbounds begin
    while(it < maxit)
      if(it>0) j=k end
    #   U[:,j] = matrixvector.stringvector(n_d,S,V[:,j],dft_s,idft_s)
      U[:,j] = H*V[:,j]
      mprod+=1
      if(it>0)
        U[:,j] = orthog(U[:,1:j-1],U[:,j])
      end
       u = norm(U[:,j])
      uinv = invcheck(u)
      U[:,j] = uinv*U[:,j]
      # Lanczos process
          while(j<m_b)
        #    v = matrixvector.adstringvector(n_d,S,U[:,j],dft_s,idft_s)
           v = adjoint(H)*U[:,j]
            mprod+=1
            v = v - u*V[:,j]
            v = orthog(V[:,1:j],v)
            fn = norm(v)
            fninv= invcheck(fn)
            v = fninv * v
                if(j<m_b-1)
                  V[:,j+1] = v
                  B[j,j] = u
                  B[j,j+1] = fn
                #  U[:,j+1] =  matrixvector.stringvector(n_d,S,V[:,j+1],dft_s,idft_s)
                 U[:,j+1] =  H*V[:,j+1]
                  mprod+=1
                  U[:,j+1] = U[:,j+1] - fn*U[:,j]
                  # ...with full reorthogonalization
                  U[:,j+1] = orthog(U[:,1:j],U[:,j+1])
                  u = norm(U[:,j+1])
                  uinv = invcheck(u)
                  U[:,j+1] = uinv * U[:,j+1]
                else
                  B[j,j] = u
   
               end
          j+=1
          end
          # end of one iteration of lanczos
         F = svd(B)
         R = fn*F.U[m_b-1,:] # Residuals
         error_end=R[nu-1]
         if(it<1)
           smax = F.S[1]
         else
           smax = max(F.S[1],smax)
         end
         # conv = count( <(tol*smax),abs.(R[1:nu]))
         conv = count( <(tol),abs.(R[1:nu]))
         # conv = k'
         if (conv<nu)
           # k = max(conv+1, k)
           k = max(conv+nu, k)
           k = min(k, m_b-3)
           k =max(k,1)
         else
           break
         end
         # resart with this part
         V[:,1:k-1] = V[:,1:m_b]*F.V[:,1:k-1]
         V[:,k] = v
         B[1:k-1,1:k-1]=Diagonal(F.S[1:k-1])
         B[1:k,k]=R[1:k]
         U[:,1:k]=U[:,1:m_b]*F.U[:,1:k]
         it+=1
       end #(while end it )
   if (it<1)
    F = svd(B)#(first svd B)
   end
    U = U[:,1:m_b]*F.U[:,1:nu]
    V = V[:,1:m_b]*F.V[:,1:nu]
   
    # println(["#mprod"              string(mprod)
    #            "#it"              string(it)
    #            "#error"              string(error_end)
    #            ])
    return (U,F.S[1:nu],V,it, mprod)
   
   end # for inbounds
   
end


# # # # # # # # # # # # # # # # # # #   truncated svd for multilevel Hankel # # # # # # #

# # # # # # # # # # # # # # # # # # #   project low-rank to Hankel # # # # # # #
function projectlowtHankel(U::AbstractArray{T},S::AbstractArray,V::AbstractArray{T},n::Array{Int64},p::Array{Int64}, q=n.-p.+1) where T
    a=ones(Int64,Tuple(p))
    b=ones(Int64,Tuple(q))
    ind_s= DSP.conv(a,b)
    Sigma=Diagonal(S)
    U=U*Sigma
    # k=size(n,1)
    r=size(U,2)
    V=conj(V)
    result=zeros(T,Tuple(n))
    for i in 1:r 
        u=reshape(U[:,i],Tuple(p))
        v=reshape(V[:,i],Tuple(q))
        result += DSP.conv(u,v)
    end
   return Hankel(result./ind_s,p)
end

# # # # # # # # # # # # # # # # # # #   project low-rank to Hankel # # # # # # #









Random.seed!(1234)
# s=rand(Int64,(100,700))
s=rand([1,2,3,4,5,6],(2000))



# s=Float64.(s)



# @time h=mulHankel(s,[50,50])


tt=Hankel(s, [1000])




# display("text/plain",h2)


x=randn(Float64,size(tt,2))
y=randn(ComplexF64,size(tt,1))



# display("text/plain",s)
# display("text/plain",x)

#  @time r1=tt*x

 
#  @time r1=tt'*y
#  tt2=tt+tt1
@time h=fullHankel(tt)
# h=Float64.(h)

# @time r2=h*x


 
#  display("text/plain", tt)
#  display("text/plain",norm(r1-r2))

 @time U_o, S_o, V_o =svd(h)
 @time U ,S, V,it,mprod=irlb(tt,5,1e-10,100,10)
#  @time U ,S, V=rimansvd(U,V,S,tt.n_d,5)
#  @time U ,S, V=riemannsvd(U,V,tt)
 @time U1 ,S1, V1,it,mprod=irlblr(tt,5,1e-10,100,10,U*0,V)
println(norm(S-S1))
oh=projectlowtHankel(U_o,S_o,V_o,tt.n_d,tt.p_d)
#  show(stdout, "text/plain", (S - S_o[1:5]))


display("text/plain", norm((tt-oh).element))
# end
