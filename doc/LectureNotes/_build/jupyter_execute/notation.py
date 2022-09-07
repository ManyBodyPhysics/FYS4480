#!/usr/bin/env python
# coding: utf-8

# <!-- HTML file automatically generated from DocOnce source (https://github.com/doconce/doconce/)
# doconce format html notation.do.txt  -->

# # Introduction to many-body physics

# ## Introduction to the course
# 
# These lecture notes aim at giving an introduction to the quantum
# mechanics of many-body systems and the methods relevant for many-body
# problems in such diverse areas as atomic, molecular, solid-state and
# nuclear physics, chemistry and materials science. A theoretical
# understanding of the behavior of quantum-mechanical many-body systems,
# that is, systems containing many interacting particles - is a
# considerable challenge in that, normally, no exact solution can be
# found.  Instead, reliable methods are needed for approximate but
# accurate simulations of such systems.
# 
# These notes cover central many-body methods used in a variety of
# fields in physics. Starting with basic definitions and second
# quantization and Wick's theorems for systems of fermions and bosons,
# we move on to the introduction of Feynman diagrams and many-body
# methods like full configuration interaction theory and various mean
# field approaches like Hartree-Fock theory and density functional
# theory. Thereafter, so-called post Hartree-Fock methods like many-body
# perturbation theory, coupled cluster theory (standard and unitary) and
# other methods will be discussed.  Finally, algorithms from quantum
# computing for solving quantum mechanical many-body problems will be
# discussed. The present set of notes contain more material than covered
# by a regular course and some of these notes can be used for
# self-studies or advanced many-body topics. Depending on the interest
# of the participants, selected methods can be emphasized over other
# ones.
# 
# This chapter serves the aim of reminding the reader of some basic
# properties of quantum mechanics, in addition to providing a brief
# review of linear algebra, which lays the foundations for most of our
# mathematical formalism. It plays a central role in our studies and
# together with computational methods and various linear algebra
# libraries, it provides us with the tools to study complicated systems
# of many interacting particles. Our short reminder on linear algebra
# allows us to catch more than one bird with a stone, we take the
# liberty of introducing some basic notations and definitions used
# throughout these lectures.

# ## Important Matrix and vector handling packages
# 
# Linear algebra plays a central role in many-body physics, as it does
# in many other mathematical descriptions of problems in physics and in
# Science in general.  There are several central software packages for
# linear algebra and eigenvalue problems. Several of the more popular
# ones have been wrapped into ofter software packages like those from
# the widely used text **Numerical Recipes**. The original source codes in
# many of the available packages are often taken from the widely used
# software package LAPACK, which follows two other popular packages
# developed in the 1970s, namely EISPACK and LINPACK.  We describe them
# shortly here.
# 
#   * LINPACK: package for linear equations and least square problems.
# 
#   * LAPACK:package for solving symmetric, unsymmetric and generalized eigenvalue problems. From LAPACK's website <http://www.netlib.org> it is possible to download for free all source codes from this library. Both C/C++ and Fortran versions are available.
# 
#   * BLAS (I, II and III): (Basic Linear Algebra Subprograms) are routines that provide standard building blocks for performing basic vector and matrix operations. Blas I is vector operations, II vector-matrix operations and III matrix-matrix operations. Highly parallelized and efficient codes, all available for download from <http://www.netlib.org>.
# 
# When dealing with matrices and vectors a central issue is memory
# handling and allocation. If our code is written in Python the way we
# declare these objects and the way they are handled, interpreted and
# used by say a linear algebra library, requires codes that interface
# our Python program with such libraries. For Python programmers,
# **Numpy** is by now the standard Python package for numerical arrays in
# Python as well as the source of functions which act on these
# arrays. These functions span from eigenvalue solvers to functions that
# compute the mean value, variance or the covariance matrix. If you are
# not familiar with how arrays are handled in say Python or compiled
# languages like C++ and Fortran, the sections in this chapter may be
# useful. For C++ programmer, **Armadillo** is widely used library for
# linear algebra and eigenvalue problems. In addition it offers a
# convenient way to handle and organize arrays. We discuss this library
# as well.  Before we proceed we believe it may be convenient to repeat
# some basic features of vectors, matrices, tensors etc. For systems of
# fermions, our ansatz for encoding the needed antisymmetry of the wave
# function, is to assume that we can build a many-body state based on
# physically relevant single-particle states. This leads the so-called
# Salter determinant.

# ## Basic Matrix Features
# 
# Matrix properties reminder

# $$
# \mathbf{A} =
#       \begin{bmatrix} a_{11} & a_{12} & a_{13} & a_{14} \\
#                                  a_{21} & a_{22} & a_{23} & a_{24} \\
#                                    a_{31} & a_{32} & a_{33} & a_{34} \\
#                                   a_{41} & a_{42} & a_{43} & a_{44}
#              \end{bmatrix}\qquad
# \mathbf{I} =
#       \begin{bmatrix} 1 & 0 & 0 & 0 \\
#                                  0 & 1 & 0 & 0 \\
#                                  0 & 0 & 1 & 0 \\
#                                  0 & 0 & 0 & 1
#              \end{bmatrix}
# $$

# The inverse of a matrix is defined by

# $$
# \mathbf{A}^{-1} \cdot \mathbf{A} = I
# $$

# <table class="dotable" border="1">
# <thead>
# <tr><th align="center">              Relations               </th> <th align="center">      Name     </th> <th align="center">                            matrix elements                            </th> </tr>
# </thead>
# <tbody>
# <tr><td align="center">   $A = A^{T}$                               </td> <td align="center">   symmetric          </td> <td align="center">   $a_{ij} = a_{ji}$                                                          </td> </tr>
# <tr><td align="center">   $A = \left (A^{T} \right )^{-1}$          </td> <td align="center">   real orthogonal    </td> <td align="center">   $\sum_k a_{ik} a_{jk} = \sum_k a_{ki} a_{kj} = \delta_{ij}$                </td> </tr>
# <tr><td align="center">   $A = A^{ * }$                             </td> <td align="center">   real matrix        </td> <td align="center">   $a_{ij} = a_{ij}^{ * }$                                                    </td> </tr>
# <tr><td align="center">   $A = A^{\dagger}$                         </td> <td align="center">   hermitian          </td> <td align="center">   $a_{ij} = a_{ji}^{ * }$                                                    </td> </tr>
# <tr><td align="center">   $A = \left (A^{\dagger} \right )^{-1}$    </td> <td align="center">   unitary            </td> <td align="center">   $\sum_k a_{ik} a_{jk}^{ * } = \sum_k a_{ki}^{ * } a_{kj} = \delta_{ij}$    </td> </tr>
# </tbody>
# </table>
# 
# Some famous Matrices we often encounter in many-body physics are
# 
#   * Diagonal if $a_{ij}=0$ for $i\ne j$
# 
#   * Upper triangular if $a_{ij}=0$ for $i > j$
# 
#   * Lower triangular if $a_{ij}=0$ for $i < j$
# 
#   * Upper Hessenberg if $a_{ij}=0$ for $i > j+1$
# 
#   * Lower Hessenberg if $a_{ij}=0$ for $i < j+1$
# 
#   * Tridiagonal if $a_{ij}=0$ for $|i -j| > 1$
# 
#   * Lower banded with bandwidth $p$: $a_{ij}=0$ for $i > j+p$
# 
#   * Upper banded with bandwidth $p$: $a_{ij}=0$ for $i < j+p$
# 
#   * Banded, block upper triangular, block lower triangular....
# 
# Furthermore, for an $N\times N$ matrix  $\mathbf{A}$ the following properties are all equivalent
# 
#   * If the inverse of $\mathbf{A}$ exists, $\mathbf{A}$ is nonsingular.
# 
#   * The equation $\mathbf{Ax}=0$ implies $\mathbf{x}=0$.
# 
#   * The rows of $\mathbf{A}$ form a basis of $R^N$.
# 
#   * The columns of $\mathbf{A}$ form a basis of $R^N$.
# 
#   * $\mathbf{A}$ is a product of elementary matrices.
# 
#   * $0$ is not eigenvalue of $\mathbf{A}$.

# ## Numpy and arrays
# 
# [Numpy](http://www.numpy.org/) provides an easy way to handle arrays in Python. The standard way to import this library is as

# In[1]:


import numpy as np
n = 10
x = np.random.normal(size=n)
print(x)


# Here we have defined a vector $x$ with $n=10$ elements with its values given by the Normal distribution $N(0,1)$.
# Another alternative is to declare a vector as follows

# In[2]:


import numpy as np
x = np.array([1, 2, 3])
print(x)


# Here we have defined a vector with three elements, with $x_0=1$, $x_1=2$ and $x_2=3$. Note that both Python and C++
# start numbering array elements from $0$ and on. This means that a vector with $n$ elements has a sequence of entities $x_0, x_1, x_2, \dots, x_{n-1}$. We could also let (recommended) Numpy to compute the logarithms of a specific array as

# In[3]:


import numpy as np
x = np.log(np.array([4, 7, 8]))
print(x)


# Here we have used Numpy's unary function $np.log$. This function is
# highly tuned to compute array elements since the code is vectorized
# and does not require looping. We normaly recommend that you use the
# Numpy intrinsic functions instead of the corresponding **log** function
# from Python's **math** module. The looping is done explicitely by the
# **np.log** function. The alternative, and slower way to compute the
# logarithms of a vector would be to write

# In[4]:


import numpy as np
from math import log
x = np.array([4, 7, 8])
for i in range(0, len(x)):
    x[i] = log(x[i])
print(x)


# We note that our code is much longer already and we need to import the **log** function from the **math** module. 
# The attentive reader will also notice that the output is $[1, 1, 2]$. Python interprets automacally our numbers as integers (like the **automatic** keyword in C++). To change this we could define our array elements to be double precision numbers as

# In[5]:


import numpy as np
x = np.log(np.array([4, 7, 8], dtype = np.float64))
print(x)


# or simply write them as double precision numbers (Python uses 64 bits as default for floating point type variables), that is

# In[6]:


import numpy as np
x = np.log(np.array([4.0, 7.0, 8.0]))
print(x)


# To check the number of bytes (remember that one byte contains eight bits for double precision variables), you can use simple use the **itemsize** functionality (the array $x$ is actually an object which inherits the functionalities defined in Numpy) as

# In[7]:


import numpy as np
x = np.log(np.array([4.0, 7.0, 8.0]))
print(x.itemsize)


# Having defined vectors, we are now ready to try out matrices. We can define a $3 \times 3 $ real matrix $\hat{A}$
# as (recall that we user lowercase letters for vectors and uppercase letters for matrices)

# In[8]:


import numpy as np
A = np.log(np.array([ [4.0, 7.0, 8.0], [3.0, 10.0, 11.0], [4.0, 5.0, 7.0] ]))
print(A)


# If we use the **shape** function we would get $(3, 3)$ as output, that is verifying that our matrix is a $3\times 3$ matrix. We can slice the matrix and print for example the first column (Python organized matrix elements in a row-major order, see below) as

# In[9]:


import numpy as np
A = np.log(np.array([ [4.0, 7.0, 8.0], [3.0, 10.0, 11.0], [4.0, 5.0, 7.0] ]))
# print the first column, row-major order and elements start with 0
print(A[:,0])


# We can continue this was by printing out other columns or rows. The example here prints out the second column

# In[10]:


import numpy as np
A = np.log(np.array([ [4.0, 7.0, 8.0], [3.0, 10.0, 11.0], [4.0, 5.0, 7.0] ]))
# print the first column, row-major order and elements start with 0
print(A[1,:])


# Numpy contains many other functionalities that allow us to slice, subdivide etc etc arrays. We strongly recommend that you look up the [Numpy website for more details](http://www.numpy.org/). Useful functions when defining a matrix are the **np.zeros** function which declares a matrix of a given dimension and sets all elements to zero

# In[11]:


import numpy as np
n = 10
# define a matrix of dimension 10 x 10 and set all elements to zero
A = np.zeros( (n, n) )
print(A)


# or initializing all elements to

# In[12]:


import numpy as np
n = 10
# define a matrix of dimension 10 x 10 and set all elements to one
A = np.ones( (n, n) )
print(A)


# or as unitarily distributed random numbers (see the material on random number generators in the statistics part)

# In[13]:


import numpy as np
n = 10
# define a matrix of dimension 10 x 10 and set all elements to random numbers with x \in [0, 1]
A = np.random.rand(n, n)
print(A)


# ## Other Matrix and Vector Operations
# 
# The following examples show how to compute various quantities like the
# **mean** value of a matrix or a vector and how to use functions like
# **reshape** and **ravel**. These are all useful quantities when scaling
# the data and preparing the data for various machine learning
# algorithms and when calculating quantities like the mean squared error
# or the variance.

# In[14]:


"""
Simple code that tests various numpy functions
"""

import numpy as np
# Simple test-matrix of dim 3 x 4
a = np.array([ [1, 2, 3], [4, 5, 6], [7, 8, 9],[10, 11, 12]],dtype=np.float64)
print(f"The test matrix:{a}")
# This is the total mean summed over all elements, which here has to be 6.5
print(f"This is the total mean summed over all elements:{np.mean(a,dtype=np.float64)}")
# This is the mean for each column, it returns an array with the mean values for each column. It returns a row-like vector
print(f"This is the mean for each column:{np.mean(a, axis=0, keepdims=True,dtype=np.float64)}")
# This is the mean value for each row, it returns an array via the keepdims option which is a column-like vector if
# keepdims=True. Else it return a row-like vector
# Try setting keepdims=False
print(f"This is the mean value for each row:{np.mean(a, axis=1, keepdims=True,dtype=np.float64)}")
# We print then the mean value for each row by  setting keepdims=False
print(f"This is the mean value for each row with keepdims false:{np.mean(a, axis=1, keepdims=False,dtype=np.float64)}")

# Ravel return a contiguous flattened array.
print(f"Flatten  the matrix:{np.ravel(a)}")
# It is the same as reshaping the matrix into a one-dimensional array
print(f"Reshape the matrix to a one-dim array:{a.reshape(-1)}")
#  ‘C’ means to index the elements in row-major, C-style order, with the last axis index changing fastest, back to the first axis index changing slowest.
# ‘F’ means to index the elements in column-major, Fortran-style order, with the first index changing fastest, and the last index changing slowest 
print(np.ravel(a, order='F'))
# When order is ‘A’, it will preserve the array’s ‘C’ or ‘F’ ordering
# ‘A’ means to read the elements in Fortran-like index order if a is Fortran contiguous in memory, C-like order otherwise.
# ‘K’ means to read the elements in the order they occur in memory, except for reversing the data when strides are negative. By default, ‘C’ index order is used.
# Transposing it
print(np.ravel(a.T))
print(np.ravel(a.T, order='A'))


# ## Definitions and vector and matrix operations
# 
# We use this section to define some central quantities as well as
# presenting our notations for single-particle and many-body states. We will
# normally tend to label vectors and states with the so-called Dirac
# **bra-ket** notation.  In general we will use lower-case and boldfaced
# letters for vectors and upper-case boldfaced letters for matrices. These ways to represent vectors will however be replaced
# with the above-mentioned Dirac **bra-ket** notation.  As example, consider the vector $\boldsymbol{x}$ defined here

# $$
# \boldsymbol{x} = \vert x \rangle = \begin{bmatrix} x_0 \\ x_1 \\ x_2 \\ \dots \\ x_{n-2} \\ x_{n-1} \end{bmatrix} \in \mathbb{C}^{n} ,
# $$

# and its hermitian conjugate

# $$
# \langle x \vert = \vert x \rangle^{\dagger}=\begin{bmatrix} x_0^* & x_1^* &  x_2^* & \dots & x_{n-2}^* & x_{n-1}^* \end{bmatrix}.
# $$

# The norm of two vectors is defined as (assuming they are not orthogonal

# $$
# \langle x \vert y \rangle =\sum_ix_i^*y_i,
# $$

# and for an orthonormal (normalized and orthogonal) basis we have

# $$
# \langle x_j \vert x_i \rangle =\delta_{ij}.
# $$

# We can also use this orthonormal basis to define the completeness relation as

# $$
# \boldsymbol{I}=\sum_{i=0}^{\infty}\vert x_i \rangle \langle _i \vert,
# $$

# where we have defined the outer product

# $$
# \vert x \rangle \langle y \vert  = \begin{bmatrix} x_0y_0^* & x_0y_1^* &  x_0y_2^* & \dots & x_0y_{n-2}^* & x_0y_{n-1}^* \\
#                                                    x_1y_0^* & x_1y_1^* &  x_1y_2^* & \dots & x_1y_{n-2}^* & x_1y_{n-1}^* \\
# 						   \dots    & \dots    & \dots     & \dots & \dots        & \dots \\                                           						         \dots    & \dots    & \dots     & \dots & \dots        & \dots \\
# 						   x_{n-2}y_0^* & x_{n-2}y_1^* &  x_{n-2}y_2^* & \dots & x_{n-2}y_{n-2}^* & x_{n-2}y_{n-1}^* \\
# 						   x_{n-1}y_0^* & x_{n-1}y_1^* &  x_{n-1}y_2^* & \dots & x_{n-1}y_{n-2}^* & x_{n-1}y_{n-1}^* \\
# \end{bmatrix}.
# $$

# The outer product defined here is going to play an important role in
# these studies as it can be used to define important projection
# operators. Since the basis of orthonormal states is, normally,
# infinite, the above sum runs over an infinite basis. However, for all
# practical reasons, unless we can find analytical solutions to a
# quantum mechanical problem, the above sum needs to be truncated.
# 
# An orthonormal basis can in turn be used to define a new basis set $\vert \phi_i \rangle$

# $$
# \vert \phi_i \rangle \sum_j C_{ij} \vert x_j\rangle,
# $$

# where the coefficients $C_{ij}=\langle \phi_i\vert x_j\rangle$ are the
# matrix elements of a transformation from a basis set $\vert \boldsymbol{x}
# \rangle$ to a new basis set $\vert \boldsymbol{\phi}\rangle$. The matrix
# $\boldsymbol{C}$ is normally a unitary matrix.
# 
# With the outer product we can define projection operators $\hat{P}$ and $\hat{Q}$.
# Consider for example the two orthonormal states

# $$
# \vert x_0 \rangle = \begin{bmatrix} 1 \\ 0 \end{bmatrix},
# $$

# and

# $$
# \vert x_1 \rangle = \begin{bmatrix} 0 \\ 1 \end{bmatrix}.
# $$

# We can easily define a new basis which is a linear combination of the two

# $$
# \vert \phi \rangle = \alpha \vert x_0 \rangle +\beta \vert x_1 \rangle
# $$

# ### Tensor products
# 
# Consider now two vectors with length $n=2$, with elements

# $$
# \vert x \rangle = \begin{bmatrix} x_0 \\ x_1 \end{bmatrix},
# $$

# and

# $$
# \vert x \rangle = \begin{bmatrix} y_0 \\ y_1 \end{bmatrix},
# $$

# $$
# \vert x \rangle \otimes \vert y \rangle = \vert xy \rangle  = \begin{bmatrix} x_0y_0 \\ x_0y_1 \\ x_1y_0 \\ x_1y_1 \end{bmatrix},
# $$

# which is now a vector of length $4$. 
# The tensor product of two matrices is defined as

# ### Representation of states
# 
# Before we proceed we need several definitions.  Throughout these
# lectures we will assume that the interacting part of the Hamiltonian
# can be approximated by a two-body interaction.  This means that our
# Hamiltonian can be written as the sum of a onebody part, which
# includes kinetic energy and an eventual external field, and a twobody
# interaction

# <!-- Equation labels as ordinary links -->
# <div id="Hnuclei"></div>
# 
# $$
# \begin{equation}
#     \hat{H} = \hat{H}_0 + \hat{H}_I 
#     = \sum_{i=1}^A \hat{h}_0(x_i) + \sum_{i < j}^A \hat{v}(r_{ij}),
# \label{Hnuclei} \tag{1}
# \end{equation}
# $$

# with

# <!-- Equation labels as ordinary links -->
# <div id="hinuclei"></div>
# 
# $$
# \begin{equation}
#   H_0=\sum_{i=1}^A \hat{h}_0(x_i).
# \label{hinuclei} \tag{2}
# \end{equation}
# $$

# The onebody part $u_{\mathrm{ext}}(x_i)$ is normally approximated by a
# harmonic oscillator potential or the Coulomb interaction an electron
# feels from the nucleus. However, other potentials are fully possible,
# such as one derived from the self-consistent solution of the
# Hartree-Fock equations to be discussed here.
# 
# Our Hamiltonian is invariant under the permutation (interchange) of
# two particles.  Since we deal with fermions however, the total wave
# function is antisymmetric.  Let $\hat{P}$ be an operator which
# interchanges two particles.  Due to the symmetries we have ascribed to
# our Hamiltonian, this operator commutes with the total Hamiltonian,

# $$
# [\hat{H},\hat{P}] = 0,
# $$

# meaning that $\Psi_{\lambda}(x_1, x_2, \dots , x_A)$ is an eigenfunction of 
# $\hat{P}$ as well, that is

# $$
# \hat{P}_{ij}\Psi_{\lambda}(x_1, x_2, \dots,x_i,\dots,x_j,\dots,x_A)=
# \beta\Psi_{\lambda}(x_1, x_2, \dots,x_i,\dots,x_j,\dots,x_A),
# $$

# where $\beta$ is the eigenvalue of $\hat{P}$. We have introduced the suffix $ij$ in order to indicate that we permute particles $i$ and $j$.
# The Pauli principle tells us that the total wave function for a system of fermions
# has to be antisymmetric, resulting in the eigenvalue $\beta = -1$.   
# 
# In our case we assume that  we can approximate the exact eigenfunction with a Slater determinant

# <!-- Equation labels as ordinary links -->
# <div id="eq:HartreeFockDet"></div>
# 
# $$
# \begin{equation}
#    \Phi(x_1, x_2,\dots ,x_A,\alpha,\beta,\dots, \sigma)=\frac{1}{\sqrt{A!}}
# \left| \begin{array}{ccccc} \psi_{\alpha}(x_1)& \psi_{\alpha}(x_2)& \dots & \dots & \psi_{\alpha}(x_A)\\
#                             \psi_{\beta}(x_1)&\psi_{\beta}(x_2)& \dots & \dots & \psi_{\beta}(x_A)\\  
#                             \dots & \dots & \dots & \dots & \dots \\
#                             \dots & \dots & \dots & \dots & \dots \\
#                      \psi_{\sigma}(x_1)&\psi_{\sigma}(x_2)& \dots & \dots & \psi_{\sigma}(x_A)\end{array} \right|, \label{eq:HartreeFockDet} \tag{3}
# \end{equation}
# $$

# where  $x_i$  stand for the coordinates and spin values of a particle $i$ and $\alpha,\beta,\dots, \gamma$ 
# are quantum numbers needed to describe remaining quantum numbers.  
# 
# Since we will deal with Fermions (identical and indistinguishable particles) we will 
# form an ansatz for a given state in terms of so-called Slater determinants determined
# by a chosen basis of single-particle functions. 
# 
# For a given $n\times n$ matrix $\mathbf{A}$ we can write its determinant

# $$
# det(\mathbf{A})=|\mathbf{A}|=
# \left| \begin{array}{ccccc} a_{11}& a_{12}& \dots & \dots & a_{1n}\\
#                             a_{21}&a_{22}& \dots & \dots & a_{2n}\\  
#                             \dots & \dots & \dots & \dots & \dots \\
#                             \dots & \dots & \dots & \dots & \dots \\
#                             a_{n1}& a_{n2}& \dots & \dots & a_{nn}\end{array} \right|,
# $$

# in a more compact form as

# $$
# |\mathbf{A}|= \sum_{i=1}^{n!}(-1)^{p_i}\hat{P}_i a_{11}a_{22}\dots a_{nn},
# $$

# where $\hat{P}_i$ is a permutation operator which permutes the column indices $1,2,3,\dots,n$
# and the sum runs over all $n!$ permutations.  The quantity $p_i$ represents the number of transpositions of column indices that are needed in order to bring a given permutation back to its initial ordering, in our case given by $a_{11}a_{22}\dots a_{nn}$ here.
# 
# A simple $2\times 2$ determinant illustrates this. We have

# $$
# det(\mathbf{A})=
# \left| \begin{array}{cc} a_{11}& a_{12}\\
#                             a_{21}&a_{22}\end{array} \right|= (-1)^0a_{11}a_{22}+(-1)^1a_{12}a_{21},
# $$

# where in the last term we have interchanged the column indices $1$ and $2$. The natural ordering we have chosen is $a_{11}a_{22}$. 
# 
# The single-particle function $\psi_{\alpha}(x_i)$  are eigenfunctions of the onebody
# Hamiltonian $h_i$, that is

# $$
# \hat{h}_0(x_i)=\hat{t}(x_i) + \hat{u}_{\mathrm{ext}}(x_i),
# $$

# with eigenvalues

# $$
# \hat{h}_0(x_i) \psi_{\alpha}(x_i)=\left(\hat{t}(x_i) + \hat{u}_{\mathrm{ext}}(x_i)\right)\psi_{\alpha}(x_i)=\varepsilon_{\alpha}\psi_{\alpha}(x_i).
# $$

# The energies $\varepsilon_{\alpha}$ are the so-called non-interacting single-particle energies, or unperturbed energies. 
# The total energy is in this case the sum over all  single-particle energies, if no two-body or more complicated
# many-body interactions are present.
# 
# Let us denote the ground state energy by $E_0$. According to the
# variational principle we have

# $$
# E_0 \le E[\Phi] = \int \Phi^*\hat{H}\Phi d\mathbf{\tau}
# $$

# where $\Phi$ is a trial function which we assume to be normalized

# $$
# \int \Phi^*\Phi d\mathbf{\tau} = 1,
# $$

# where we have used the shorthand $d\mathbf{\tau}=dx_1dr_2\dots dr_A$.
# 
# In the Hartree-Fock method the trial function is the Slater
# determinant of Eq. ([3](#eq:HartreeFockDet)) which can be rewritten as

# $$
# \Phi(x_1,x_2,\dots,x_A,\alpha,\beta,\dots,\nu) = \frac{1}{\sqrt{A!}}\sum_{P} (-)^P\hat{P}\psi_{\alpha}(x_1)
#     \psi_{\beta}(x_2)\dots\psi_{\nu}(x_A)=\sqrt{A!}\hat{A}\Phi_H,
# $$

# where we have introduced the antisymmetrization operator $\hat{A}$ defined by the 
# summation over all possible permutations of two particles.
# 
# It is defined as

# <!-- Equation labels as ordinary links -->
# <div id="antiSymmetryOperator"></div>
# 
# $$
# \begin{equation}
#   \hat{A} = \frac{1}{A!}\sum_{p} (-)^p\hat{P},
# \label{antiSymmetryOperator} \tag{4}
# \end{equation}
# $$

# with $p$ standing for the number of permutations. We have introduced for later use the so-called
# Hartree-function, defined by the simple product of all possible single-particle functions

# $$
# \Phi_H(x_1,x_2,\dots,x_A,\alpha,\beta,\dots,\nu) =
#   \psi_{\alpha}(x_1)
#     \psi_{\beta}(x_2)\dots\psi_{\nu}(x_A).
# $$

# Both $\hat{H}_0$ and $\hat{H}_I$ are invariant under all possible permutations of any two particles
# and hence commute with $\hat{A}$

# <!-- Equation labels as ordinary links -->
# <div id="commutionAntiSym"></div>
# 
# $$
# \begin{equation}
#   [H_0,\hat{A}] = [H_I,\hat{A}] = 0. \label{commutionAntiSym} \tag{5}
# \end{equation}
# $$

# Furthermore, $\hat{A}$ satisfies

# <!-- Equation labels as ordinary links -->
# <div id="AntiSymSquared"></div>
# 
# $$
# \begin{equation}
#   \hat{A}^2 = \hat{A},  \label{AntiSymSquared} \tag{6}
# \end{equation}
# $$

# since every permutation of the Slater
# determinant reproduces it. 
# 
# The expectation value of $\hat{H}_0$

# $$
# \int \Phi^*\hat{H}_0\Phi d\mathbf{\tau} 
#   = A! \int \Phi_H^*\hat{A}\hat{H}_0\hat{A}\Phi_H d\mathbf{\tau}
# $$

# is readily reduced to

# $$
# \int \Phi^*\hat{H}_0\Phi d\mathbf{\tau} 
#   = A! \int \Phi_H^*\hat{H}_0\hat{A}\Phi_H d\mathbf{\tau},
# $$

# where we have used Eqs. ([5](#commutionAntiSym)) and
# ([6](#AntiSymSquared)). The next step is to replace the antisymmetrization
# operator by its definition and to
# replace $\hat{H}_0$ with the sum of one-body operators

# $$
# \int \Phi^*\hat{H}_0\Phi  d\mathbf{\tau}
#   = \sum_{i=1}^A \sum_{p} (-)^p\int 
#   \Phi_H^*\hat{h}_0\hat{P}\Phi_H d\mathbf{\tau}.
# $$

# The integral vanishes if two or more particles are permuted in only one
# of the Hartree-functions $\Phi_H$ because the individual single-particle wave functions are
# orthogonal. We obtain then

# $$
# \int \Phi^*\hat{H}_0\Phi  d\mathbf{\tau}= \sum_{i=1}^A \int \Phi_H^*\hat{h}_0\Phi_H  d\mathbf{\tau}.
# $$

# Orthogonality of the single-particle functions allows us to further simplify the integral, and we
# arrive at the following expression for the expectation values of the
# sum of one-body Hamiltonians

# <!-- Equation labels as ordinary links -->
# <div id="H1Expectation"></div>
# 
# $$
# \begin{equation}
#   \int \Phi^*\hat{H}_0\Phi  d\mathbf{\tau}
#   = \sum_{\mu=1}^A \int \psi_{\mu}^*(x)\hat{h}_0\psi_{\mu}(x)dx
#   d\mathbf{r}.
# \label{H1Expectation} \tag{7}
# \end{equation}
# $$

# We introduce the following shorthand for the above integral

# $$
# \langle \mu | \hat{h}_0 | \mu \rangle = \int \psi_{\mu}^*(x)\hat{h}_0\psi_{\mu}(x)dx,
# $$

# and rewrite Eq. ([7](#H1Expectation)) as

# <!-- Equation labels as ordinary links -->
# <div id="H1Expectation1"></div>
# 
# $$
# \begin{equation}
#   \int \Phi^*\hat{H}_0\Phi  d\tau
#   = \sum_{\mu=1}^A \langle \mu | \hat{h}_0 | \mu \rangle.
# \label{H1Expectation1} \tag{8}
# \end{equation}
# $$

# The expectation value of the two-body part of the Hamiltonian is obtained in a
# similar manner. We have

# $$
# \int \Phi^*\hat{H}_I\Phi d\mathbf{\tau} 
#   = A! \int \Phi_H^*\hat{A}\hat{H}_I\hat{A}\Phi_H d\mathbf{\tau},
# $$

# which reduces to

# $$
# \int \Phi^*\hat{H}_I\Phi d\mathbf{\tau} 
#   = \sum_{i\le j=1}^A \sum_{p} (-)^p\int 
#   \Phi_H^*\hat{v}(r_{ij})\hat{P}\Phi_H d\mathbf{\tau},
# $$

# by following the same arguments as for the one-body
# Hamiltonian. 
# 
# Because of the dependence on the inter-particle distance $r_{ij}$,  permutations of
# any two particles no longer vanish, and we get

# $$
# \int \Phi^*\hat{H}_I\Phi d\mathbf{\tau} 
#   = \sum_{i < j=1}^A \int  
#   \Phi_H^*\hat{v}(r_{ij})(1-P_{ij})\Phi_H d\mathbf{\tau}.
# $$

# where $P_{ij}$ is the permutation operator that interchanges
# particle $i$ and particle $j$. Again we use the assumption that the single-particle wave functions
# are orthogonal. 
# 
# We obtain

# <!-- Equation labels as ordinary links -->
# <div id="_auto1"></div>
# 
# $$
# \begin{equation}
#   \int \Phi^*\hat{H}_I\Phi d\mathbf{\tau} 
#   = \frac{1}{2}\sum_{\mu=1}^A\sum_{\nu=1}^A
#     \left[ \int \psi_{\mu}^*(x_i)\psi_{\nu}^*(x_j)\hat{v}(r_{ij})\psi_{\mu}(x_i)\psi_{\nu}(x_j)
#     dx_idx_j \right.
# \label{_auto1} \tag{9}
# \end{equation}
# $$

# <!-- Equation labels as ordinary links -->
# <div id="H2Expectation"></div>
# 
# $$
# \begin{equation} 
#   \left.
#   - \int \psi_{\mu}^*(x_i)\psi_{\nu}^*(x_j)
#   \hat{v}(r_{ij})\psi_{\nu}(x_i)\psi_{\mu}(x_j)
#   dx_idx_j
#   \right]. \label{H2Expectation} \tag{10}
# \end{equation}
# $$

# The first term is the so-called direct term. It is frequently also called the  Hartree term, 
# while the second is due to the Pauli principle and is called
# the exchange term or just the Fock term.
# The factor  $1/2$ is introduced because we now run over
# all pairs twice. 
# 
# The last equation allows us to  introduce some further definitions.  
# The single-particle wave functions $\psi_{\mu}(x)$, defined by the quantum numbers $\mu$ and $x$
# are defined as the overlap

# $$
# \psi_{\alpha}(x)  = \langle x | \alpha \rangle .
# $$

# We introduce the following shorthands for the above two integrals

# $$
# \langle \mu\nu|\hat{v}|\mu\nu\rangle =  \int \psi_{\mu}^*(x_i)\psi_{\nu}^*(x_j)\hat{v}(r_{ij})\psi_{\mu}(x_i)\psi_{\nu}(x_j)
#     dx_idx_j,
# $$

# and

# $$
# \langle \mu\nu|\hat{v}|\nu\mu\rangle = \int \psi_{\mu}^*(x_i)\psi_{\nu}^*(x_j)
#   \hat{v}(r_{ij})\psi_{\nu}(x_i)\psi_{\mu}(x_j)
#   dx_idx_j.
# $$

# ## Preparing for later studies: varying the coefficients of a wave function expansion and orthogonal transformations
# 
# It is common to  expand the single-particle functions in a known basis  and vary the coefficients, 
# that is, the new single-particle wave function is written as a linear expansion
# in terms of a fixed chosen orthogonal basis (for example the well-known harmonic oscillator functions or the hydrogen-like functions etc).
# We define our new single-particle basis (this is a normal approach for Hartree-Fock theory) by performing a unitary transformation 
# on our previous basis (labelled with greek indices) as

# <!-- Equation labels as ordinary links -->
# <div id="eq:newbasis"></div>
# 
# $$
# \begin{equation}
# \psi_p^{new}  = \sum_{\lambda} C_{p\lambda}\phi_{\lambda}. \label{eq:newbasis} \tag{11}
# \end{equation}
# $$

# In this case we vary the coefficients $C_{p\lambda}$. If the basis has infinitely many solutions, we need
# to truncate the above sum.  We assume that the basis $\phi_{\lambda}$ is orthogonal.
# 
# It is normal to choose a single-particle basis defined as the eigenfunctions
# of parts of the full Hamiltonian. The typical situation consists of the solutions of the one-body part of the Hamiltonian, that is we have

# $$
# \hat{h}_0\phi_{\lambda}=\epsilon_{\lambda}\phi_{\lambda}.
# $$

# The single-particle wave functions $\phi_{\lambda}(\mathbf{r})$, defined by the quantum numbers $\lambda$ and $\mathbf{r}$
# are defined as the overlap

# $$
# \phi_{\lambda}(\mathbf{r})  = \langle \mathbf{r} | \lambda \rangle .
# $$

# In deriving the Hartree-Fock equations, we  will expand the single-particle functions in a known basis  and vary the coefficients, 
# that is, the new single-particle wave function is written as a linear expansion
# in terms of a fixed chosen orthogonal basis (for example the well-known harmonic oscillator functions or the hydrogen-like functions etc).
# 
# We stated that a unitary transformation keeps the orthogonality. To see this consider first a basis of vectors $\mathbf{v}_i$,

# $$
# \mathbf{v}_i = \begin{bmatrix} v_{i1} \\ \dots \\ \dots \\v_{in} \end{bmatrix}
# $$

# We assume that the basis is orthogonal, that is

# $$
# \mathbf{v}_j^T\mathbf{v}_i = \delta_{ij}.
# $$

# An orthogonal or unitary transformation

# $$
# \mathbf{w}_i=\mathbf{U}\mathbf{v}_i,
# $$

# preserves the dot product and orthogonality since

# $$
# \mathbf{w}_j^T\mathbf{w}_i=(\mathbf{U}\mathbf{v}_j)^T\mathbf{U}\mathbf{v}_i=\mathbf{v}_j^T\mathbf{U}^T\mathbf{U}\mathbf{v}_i= \mathbf{v}_j^T\mathbf{v}_i = \delta_{ij}.
# $$

# This means that if the coefficients $C_{p\lambda}$ belong to a unitary or orthogonal trasformation (using the Dirac bra-ket notation)

# $$
# \vert p\rangle  = \sum_{\lambda} C_{p\lambda}\vert\lambda\rangle,
# $$

# orthogonality is preserved, that is $\langle \alpha \vert \beta\rangle = \delta_{\alpha\beta}$
# and $\langle p \vert q\rangle = \delta_{pq}$. 
# 
# This propertry is extremely useful when we build up a basis of many-body Stater determinant based states. 
# 
# **Note also that although a basis $\vert \alpha\rangle$ contains an infinity of states, for practical calculations we have always to make some truncations.** 
# 
# Before we develop for example the Hartree-Fock equations, there is another very useful property of determinants that we will use both in connection with Hartree-Fock calculations and later shell-model calculations.  
# 
# Consider the following determinant

# $$
# \left| \begin{array}{cc} \alpha_1b_{11}+\alpha_2sb_{12}& a_{12}\\
#                          \alpha_1b_{21}+\alpha_2b_{22}&a_{22}\end{array} \right|=\alpha_1\left|\begin{array}{cc} b_{11}& a_{12}\\
#                          b_{21}&a_{22}\end{array} \right|+\alpha_2\left| \begin{array}{cc} b_{12}& a_{12}\\b_{22}&a_{22}\end{array} \right|
# $$

# We can generalize this to  an $n\times n$ matrix and have

# $$
# \left| \begin{array}{cccccc} a_{11}& a_{12} & \dots & \sum_{k=1}^n c_k b_{1k} &\dots & a_{1n}\\
# a_{21}& a_{22} & \dots & \sum_{k=1}^n c_k b_{2k} &\dots & a_{2n}\\
# \dots & \dots & \dots & \dots & \dots & \dots \\
# \dots & \dots & \dots & \dots & \dots & \dots \\
# a_{n1}& a_{n2} & \dots & \sum_{k=1}^n c_k b_{nk} &\dots & a_{nn}\end{array} \right|=
# \sum_{k=1}^n c_k\left| \begin{array}{cccccc} a_{11}& a_{12} & \dots &  b_{1k} &\dots & a_{1n}\\
# a_{21}& a_{22} & \dots &  b_{2k} &\dots & a_{2n}\\
# \dots & \dots & \dots & \dots & \dots & \dots\\
# \dots & \dots & \dots & \dots & \dots & \dots\\
# a_{n1}& a_{n2} & \dots &  b_{nk} &\dots & a_{nn}\end{array} \right| .
# $$

# This is a property we will use in our Hartree-Fock discussions. 
# 
# We can generalize the previous results, now 
# with all elements $a_{ij}$  being given as functions of 
# linear combinations  of various coefficients $c$ and elements $b_{ij}$,

# $$
# \left| \begin{array}{cccccc} \sum_{k=1}^n b_{1k}c_{k1}& \sum_{k=1}^n b_{1k}c_{k2} & \dots & \sum_{k=1}^n b_{1k}c_{kj}  &\dots & \sum_{k=1}^n b_{1k}c_{kn}\\
# \sum_{k=1}^n b_{2k}c_{k1}& \sum_{k=1}^n b_{2k}c_{k2} & \dots & \sum_{k=1}^n b_{2k}c_{kj} &\dots & \sum_{k=1}^n b_{2k}c_{kn}\\
# \dots & \dots & \dots & \dots & \dots & \dots \\
# \dots & \dots & \dots & \dots & \dots &\dots \\
# \sum_{k=1}^n b_{nk}c_{k1}& \sum_{k=1}^n b_{nk}c_{k2} & \dots & \sum_{k=1}^n b_{nk}c_{kj} &\dots & \sum_{k=1}^n b_{nk}c_{kn}\end{array} \right|=det(\mathbf{C})det(\mathbf{B}),
# $$

# where $det(\mathbf{C})$ and $det(\mathbf{B})$ are the determinants of $n\times n$ matrices
# with elements $c_{ij}$ and $b_{ij}$ respectively.  
# This is a property we will use in our Hartree-Fock discussions. Convince yourself about the correctness of the above expression by setting $n=2$. 
# 
# With our definition of the new basis in terms of an orthogonal basis we have

# $$
# \psi_p(x)  = \sum_{\lambda} C_{p\lambda}\phi_{\lambda}(x).
# $$

# If the coefficients $C_{p\lambda}$ belong to an orthogonal or unitary matrix, the new basis
# is also orthogonal. 
# Our Slater determinant in the new basis $\psi_p(x)$ is written as

# $$
# \frac{1}{\sqrt{A!}}
# \left| \begin{array}{ccccc} \psi_{p}(x_1)& \psi_{p}(x_2)& \dots & \dots & \psi_{p}(x_A)\\
#                             \psi_{q}(x_1)&\psi_{q}(x_2)& \dots & \dots & \psi_{q}(x_A)\\  
#                             \dots & \dots & \dots & \dots & \dots \\
#                             \dots & \dots & \dots & \dots & \dots \\
#                      \psi_{t}(x_1)&\psi_{t}(x_2)& \dots & \dots & \psi_{t}(x_A)\end{array} \right|=\frac{1}{\sqrt{A!}}
# \left| \begin{array}{ccccc} \sum_{\lambda} C_{p\lambda}\phi_{\lambda}(x_1)& \sum_{\lambda} C_{p\lambda}\phi_{\lambda}(x_2)& \dots & \dots & \sum_{\lambda} C_{p\lambda}\phi_{\lambda}(x_A)\\
#                             \sum_{\lambda} C_{q\lambda}\phi_{\lambda}(x_1)&\sum_{\lambda} C_{q\lambda}\phi_{\lambda}(x_2)& \dots & \dots & \sum_{\lambda} C_{q\lambda}\phi_{\lambda}(x_A)\\  
#                             \dots & \dots & \dots & \dots & \dots \\
#                             \dots & \dots & \dots & \dots & \dots \\
#                      \sum_{\lambda} C_{t\lambda}\phi_{\lambda}(x_1)&\sum_{\lambda} C_{t\lambda}\phi_{\lambda}(x_2)& \dots & \dots & \sum_{\lambda} C_{t\lambda}\phi_{\lambda}(x_A)\end{array} \right|,
# $$

# which is nothing but $det(\mathbf{C})det(\Phi)$, with $det(\Phi)$ being the determinant given by the basis functions $\phi_{\lambda}(x)$. 
# 
# In our discussions hereafter we will use our definitions of single-particle states above and below the Fermi ($F$) level given by the labels
# $ijkl\dots \le F$ for so-called single-hole states and $abcd\dots > F$ for so-called particle states.
# For general single-particle states we employ the labels $pqrs\dots$. 
# 
# The energy functional is

# $$
# E[\Phi] 
#   = \sum_{\mu=1}^A \langle \mu | h | \mu \rangle +
#   \frac{1}{2}\sum_{{\mu}=1}^A\sum_{{\nu}=1}^A \langle \mu\nu|\hat{v}|\mu\nu\rangle_{AS},
# $$

# we found the expression for the energy functional in terms of the basis function $\phi_{\lambda}(\mathbf{r})$. We then  varied the above energy functional with respect to the basis functions $|\mu \rangle$. 
# Now we are interested in defining a new basis defined in terms of
# a chosen basis as defined in Eq. ([11](#eq:newbasis)). We can then rewrite the energy functional as

# <!-- Equation labels as ordinary links -->
# <div id="FunctionalEPhi2"></div>
# 
# $$
# \begin{equation}
#   E[\Phi^{New}] 
#   = \sum_{i=1}^A \langle i | h | i \rangle +
#   \frac{1}{2}\sum_{ij=1}^A\langle ij|\hat{v}|ij\rangle_{AS}, \label{FunctionalEPhi2} \tag{12}
# \end{equation}
# $$

# where $\Phi^{New}$ is the new Slater determinant defined by the new basis of Eq. ([11](#eq:newbasis)). 
# 
# Using Eq. ([11](#eq:newbasis)) we can rewrite Eq. ([12](#FunctionalEPhi2)) as

# <!-- Equation labels as ordinary links -->
# <div id="FunctionalEPhi3"></div>
# 
# $$
# \begin{equation}
#   E[\Psi] 
#   = \sum_{i=1}^A \sum_{\alpha\beta} C^*_{i\alpha}C_{i\beta}\langle \alpha | h | \beta \rangle +
#   \frac{1}{2}\sum_{ij=1}^A\sum_{{\alpha\beta\gamma\delta}} C^*_{i\alpha}C^*_{j\beta}C_{i\gamma}C_{j\delta}\langle \alpha\beta|\hat{v}|\gamma\delta\rangle_{AS}. \label{FunctionalEPhi3} \tag{13}
# \end{equation}
# $$
