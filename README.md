# defmm
defmm (Directional Equispaced interpolation-based Fast Multipole Method) is a header only C++ library. This is a light version of the code presented in https://tel.archives-ouvertes.fr/tel-03203231/

The goal of this library is to efficiently process any particle distribution with possibly oscillatory kernels by using directional interpolation.
Far field blocks are evaluated fastly by exploiting Fourier techniques.

## Dependencies
defmm uses OpenMP, BLAS and FFTW3.

## Install and test
(Using g++) In the root directory, run

$ make testinterface

$ ./testinterface

This simple code runs the FMM on a small non-uniform particle distribution in the high-frequeny regime.
Input files list the coordinates of particles using the following format:

x y z

x y z

[...]

## Licence
The defmm library is under Lesser Gnu Public Licence (LGPL).

## Authors
Authors: Igor Chollet (Sorbonne Université / Inria), Xavier Claeys (Sorbonne Université / Inria), Pierre Fortin (Sorbonne Université) and Laura Grigori (Sorbonne Université / Inria).

## Credits
We acknowledge support of Institut des Sciences du Calcul et des Données (Sorbonne Université) as well as the support of the Inria Paris (Alpines team) and of the French National Research Agency (ANR) under the grant ANR-15-CE23-0017-01.

![Image](https://github.com/xclaeys/BemTool/blob/master/doc/Logo-anr.png?raw=true)
