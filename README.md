# multiblockHankelMatrices.jl
===========

This is a Julia package related to Multiblock Hankel Matrices structures. It mainly used for spectral compressed seeing and has several functions implemented with FFT.

Fast matrix multiplication for multiblock Hankel matrices in Julia

## Note

Multiplication of large matrices for multiblock Hankel matrices
are computed with FFTs.
To be able to use these methods, you have to install and load a package that implements
the [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl) interface such
as [FFTW.jl](https://github.com/JuliaMath/FFTW.jl):

```julia
using FFTW
```

## Introduction


### Hankel
A Hankel matrix has constant anti-diagonals, for example,
```julia
[ 1  2  3
  2  3  4
  3  4  5 ]
```
We define the Hankel matrix via the array involved and the index of the pencil parameter: 
```julia
Hankel(s,p)
```
where s = [1 2 3 4 5 ] is the array to store the elements for the first column and last row, and p = 3 is the pencil parameter. 

### Multiblock  Hankel 

The reason we set the pencil parameter is easily expanded to multiblock Hankel matrices, i.e 2-D Hankel: 
```julia
[1 2 3  2 3 4
 2 3 4  3 4 5
 3 4 5  4 5 6
 2 3 4  3 4 5 
 3 4 5  4 5 6
 4 5 6  5 6 7]
```
where s=[ 1 2 3 4 5; 2 3 4 4 6; 3 4 5 6 7; ], p=[3,2]
The scalar multiplication, `conj`, `+`, `-` and `transpose` are implemented with dense form. 
Unless the fullHankel function is called, the explicit matrix form is not applied for computational and storage efficiency.
```julia
fullHankel(Hankel(s,p))
```


## functions

### Matrix-vector multiplication 
The FFT is implemented for fast matrix-vector multiplication unless the size of the one-dimensional Hankel matrix is less than 512. Moreover, the initialization for the required arrays FFT is planed in advance only once with `LinearAlgebra.factorize` to avoid repetitive calculations.

### Truncated SVD for Hankel matrix



### Adjoint operator $$ \mathcal{H}^* $$ on low-rank matrix 
