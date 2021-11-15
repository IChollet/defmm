//===================================================================
//
// Authors: Igor Chollet, Xavier Claeys, Pierre Fortin, Laura Grigori
//
//  Copyright (C) 2020 Sorbonne Universite, Inria
//  All right reserved.
//
//  This file is part of defmm.
//
//  defmm is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  defmm is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  (see ../../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with defmm.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================

#ifndef BLAS_PARSER_CONSTS
#define BLAS_PARSER_CONSTS
#define BLAS_S float
#define BLAS_D double
#define BLAS_C std::complex<float>
#define BLAS_Z std::complex<double>
BLAS_S S_ZERO     =  0.0;
BLAS_S S_ONE      =  1.0;
BLAS_S S_MONE     = -1.0;
BLAS_D D_ZERO     =  0.0;
BLAS_D D_ONE      =  1.0;
BLAS_D D_MONE     = -1.0;
BLAS_C C_ZERO     =  0.0;
BLAS_C C_ONE      =  1.0;
BLAS_C C_MONE     = -1.0;
BLAS_Z Z_ZERO     =  0.0;
BLAS_Z Z_ONE      =  1.0;
BLAS_Z Z_MONE     = -1.0;
int    IN_ONE     =  1  ;
const char* charN = "N" ;
const char* charT = "T" ;
const char* charC = "C" ;
const char* charS = "S" ;
#endif
