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

#ifndef GAIA_ALIASES_HPP
#define GAIA_ALIASES_HPP

#include "gaia_vectors.hpp"
#include "gaia_small_vectors.hpp"
#include "gaia_matrices.hpp"

// typedefs
typedef double            flt;
typedef std::complex<flt> cplx;
typedef gaia::vec<3,flt>  vec3;
typedef gaia::vec<3,int>  ivec3;
typedef gaia::Vec<flt>    Vecf;
typedef gaia::Mat<flt>    Matf;
typedef gaia::Vec<cplx>   Vecc;
typedef gaia::Mat<cplx>   Matc;
typedef gaia::Vec<int>    Veci;
typedef gaia::Mat<int>    Mati;
typedef gaia::Vec<bool>   Vecb;
typedef gaia::Mat<bool>   Matb;
template <int DIM> using veci = gaia::vec<DIM,int>;
template <int DIM> using vecf = gaia::vec<DIM,flt>;
template <int DIM> using vecc = gaia::vec<DIM,cplx>;
template <typename U> using Vec = gaia::Vec<U>;
template <typename U> using Mat = gaia::Mat<U>; 

// primitives usuelles
const cplx cplxi = cplx(0.,1.);
const flt  pi4   = 4.*M_PI;
const cplx cplx0 = cplx(0.,0.);
#define urand rand()/flt(RAND_MAX)

#endif
