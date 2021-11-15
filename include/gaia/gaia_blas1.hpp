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

#ifndef GAIA_BLAS_1_HPP
#define GAIA_BLAS_1_HPP

#include <complex>
#include "gaia_vectors.hpp"
#include "gaia_blas_consts.hpp"

namespace gaia{
  extern "C"{
    BLAS_S sdot_ (const int*, BLAS_S*, int*, BLAS_S*, int*);
    BLAS_D ddot_ (const int*, BLAS_D*, int*, BLAS_D*, int*);
    BLAS_C cdotc_(const int*, BLAS_C*, int*, BLAS_C*, int*);
    BLAS_Z zdotc_(const int*, BLAS_Z*, int*, BLAS_Z*, int*);
    BLAS_C cdotu_(const int*, BLAS_C*, int*, BLAS_C*, int*);
    BLAS_Z zdotu_(const int*, BLAS_Z*, int*, BLAS_Z*, int*);
  }
  inline BLAS_S dot (Vec<BLAS_S>& v1, Vec<BLAS_S>& v2){return sdot_ (pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_D dot (Vec<BLAS_D>& v1, Vec<BLAS_D>& v2){return ddot_ (pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_C dotc(Vec<BLAS_C>& v1, Vec<BLAS_C>& v2){return cdotc_(pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_Z dotc(Vec<BLAS_Z>& v1, Vec<BLAS_Z>& v2){return zdotc_(pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_C dotu(Vec<BLAS_C>& v1, Vec<BLAS_C>& v2){return cdotu_(pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_Z dotu(Vec<BLAS_Z>& v1, Vec<BLAS_Z>& v2){return zdotu_(pSize(v1), K_(v1), &IN_ONE, K_(v2), &IN_ONE);}
  inline BLAS_D dot (Vec<BLAS_D>& v1, Vec<BLAS_D>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return ddot_ (&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  inline BLAS_S dot (Vec<BLAS_S>& v1, Vec<BLAS_S>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return sdot_ (&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  inline BLAS_C dotc(Vec<BLAS_C>& v1, Vec<BLAS_C>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return cdotc_(&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  inline BLAS_Z dotc(Vec<BLAS_Z>& v1, Vec<BLAS_Z>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return zdotc_(&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  inline BLAS_C dotu(Vec<BLAS_C>& v1, Vec<BLAS_C>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return cdotu_(&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  inline BLAS_Z dotu(Vec<BLAS_Z>& v1, Vec<BLAS_Z>& v2, int startx, int stepx, int starty, int stepy, int nbdata){return zdotu_(&nbdata, K_(v1)+startx, &stepx, K_(v2)+starty, &stepy);}
  extern "C"{
    void saxpy_(const int*, BLAS_S*, BLAS_S*, int*, BLAS_S*, int*);
    void daxpy_(const int*, BLAS_D*, BLAS_D*, int*, BLAS_D*, int*);
    void caxpy_(const int*, BLAS_C*, BLAS_C*, int*, BLAS_C*, int*);
    void zaxpy_(const int*, BLAS_Z*, BLAS_Z*, int*, BLAS_Z*, int*);
  }
  inline void axpy(BLAS_S a, Vec<BLAS_S>& x, Vec<BLAS_S>& y){saxpy_(pSize(x), &a, K_(x), &IN_ONE, K_(y), &IN_ONE);}
  inline void axpy(BLAS_D a, Vec<BLAS_D>& x, Vec<BLAS_D>& y){daxpy_(pSize(x), &a, K_(x), &IN_ONE, K_(y), &IN_ONE);}
  inline void axpy(BLAS_C a, Vec<BLAS_C>& x, Vec<BLAS_C>& y){caxpy_(pSize(x), &a, K_(x), &IN_ONE, K_(y), &IN_ONE);}
  inline void axpy(BLAS_Z a, Vec<BLAS_Z>& x, Vec<BLAS_Z>& y){zaxpy_(pSize(x), &a, K_(x), &IN_ONE, K_(y), &IN_ONE);}
  inline void axpy(BLAS_S a, Vec<BLAS_S>& x, Vec<BLAS_S>& y, int startx, int stepx, int starty, int stepy, int nbdata){saxpy_(&nbdata, &a, K_(x)+startx, &stepx, K_(y)+starty, &stepy);}
  inline void axpy(BLAS_D a, Vec<BLAS_D>& x, Vec<BLAS_D>& y, int startx, int stepx, int starty, int stepy, int nbdata){daxpy_(&nbdata, &a, K_(x)+startx, &stepx, K_(y)+starty, &stepy);}
  inline void axpy(BLAS_C a, Vec<BLAS_C>& x, Vec<BLAS_C>& y, int startx, int stepx, int starty, int stepy, int nbdata){caxpy_(&nbdata, &a, K_(x)+startx, &stepx, K_(y)+starty, &stepy);}
  inline void axpy(BLAS_Z a, Vec<BLAS_Z>& x, Vec<BLAS_Z>& y, int startx, int stepx, int starty, int stepy, int nbdata){zaxpy_(&nbdata, &a, K_(x)+startx, &stepx, K_(y)+starty, &stepy);}
  extern "C"{
    void sscal_(const int*, BLAS_S*, BLAS_S*, int*);
    void dscal_(const int*, BLAS_D*, BLAS_D*, int*);
    void cscal_(const int*, BLAS_C*, BLAS_C*, int*);
    void zscal_(const int*, BLAS_Z*, BLAS_Z*, int*);
  }
  inline void scal(BLAS_S a, Vec<BLAS_S>& x){sscal_(pSize(x),&a,K_(x),&IN_ONE);}
  inline void scal(BLAS_D a, Vec<BLAS_D>& x){dscal_(pSize(x),&a,K_(x),&IN_ONE);}
  inline void scal(BLAS_C a, Vec<BLAS_C>& x){cscal_(pSize(x),&a,K_(x),&IN_ONE);}
  inline void scal(BLAS_Z a, Vec<BLAS_Z>& x){zscal_(pSize(x),&a,K_(x),&IN_ONE);}
  inline void scal(BLAS_S a, Vec<BLAS_S>& x, int startx, int stepx, int nbdata){sscal_(&nbdata,&a,K_(x)+startx,&stepx);}
  inline void scal(BLAS_D a, Vec<BLAS_D>& x, int startx, int stepx, int nbdata){dscal_(&nbdata,&a,K_(x)+startx,&stepx);}
  inline void scal(BLAS_C a, Vec<BLAS_C>& x, int startx, int stepx, int nbdata){cscal_(&nbdata,&a,K_(x)+startx,&stepx);}
  inline void scal(BLAS_Z a, Vec<BLAS_Z>& x, int startx, int stepx, int nbdata){zscal_(&nbdata,&a,K_(x)+startx,&stepx);}
  extern "C"{
    void scopy_(const int*, BLAS_S*, int*, BLAS_S*, int*);
    void dcopy_(const int*, BLAS_D*, int*, BLAS_D*, int*);
    void ccopy_(const int*, BLAS_C*, int*, BLAS_C*, int*);
    void zcopy_(const int*, BLAS_Z*, int*, BLAS_Z*, int*);
  }
  inline void copy(Vec<BLAS_S>& x, Vec<BLAS_S>& y){scopy_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void copy(Vec<BLAS_D>& x, Vec<BLAS_D>& y){dcopy_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void copy(Vec<BLAS_C>& x, Vec<BLAS_C>& y){ccopy_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void copy(Vec<BLAS_Z>& x, Vec<BLAS_Z>& y){zcopy_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void copy(Vec<BLAS_S>& x, Vec<BLAS_S>& y, int startx, int stepx, int starty, int stepy, int nbdata){scopy_(&nbdata,K_(x)+startx,&stepx,K_(y)+starty,&stepy);}
  inline void copy(Vec<BLAS_D>& x, Vec<BLAS_D>& y, int startx, int stepx, int starty, int stepy, int nbdata){dcopy_(&nbdata,K_(x)+startx,&stepx,K_(y)+starty,&stepy);}
  inline void copy(Vec<BLAS_C>& x, Vec<BLAS_C>& y, int startx, int stepx, int starty, int stepy, int nbdata){ccopy_(&nbdata,K_(x)+startx,&stepx,K_(y)+starty,&stepy);}
  inline void copy(Vec<BLAS_Z>& x, Vec<BLAS_Z>& y, int startx, int stepx, int starty, int stepy, int nbdata){zcopy_(&nbdata,K_(x)+startx,&stepx,K_(y)+starty,&stepy);}
  extern "C"{
    void sswap_(const int*, BLAS_S*, int*, BLAS_S*, int*);
    void dswap_(const int*, BLAS_D*, int*, BLAS_D*, int*);
    void cswap_(const int*, BLAS_C*, int*, BLAS_C*, int*);
    void zswap_(const int*, BLAS_Z*, int*, BLAS_Z*, int*);
  }
  inline void swap(Vec<BLAS_S>& x, Vec<BLAS_S>& y){sswap_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void swap(Vec<BLAS_D>& x, Vec<BLAS_D>& y){dswap_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void swap(Vec<BLAS_C>& x, Vec<BLAS_C>& y){cswap_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  inline void swap(Vec<BLAS_Z>& x, Vec<BLAS_Z>& y){zswap_(pSize(x),K_(x),&IN_ONE,K_(y),&IN_ONE);}
  extern "C"{
    BLAS_S snrm2_ (const int*, BLAS_S*, int*);
    BLAS_D dnrm2_ (const int*, BLAS_D*, int*);
    BLAS_S scnrm2_(const int*, BLAS_C*, int*);
    BLAS_D dznrm2_(const int*, BLAS_Z*, int*);
  }
  inline BLAS_S nrm2(Vec<BLAS_S>& x){return snrm2_ (pSize(x),K_(x),&IN_ONE);}
  inline BLAS_D nrm2(Vec<BLAS_D>& x){return dnrm2_ (pSize(x),K_(x),&IN_ONE);}
  inline BLAS_S nrm2(Vec<BLAS_C>& x){return scnrm2_(pSize(x),K_(x),&IN_ONE);}
  inline BLAS_D nrm2(Vec<BLAS_Z>& x){return dznrm2_(pSize(x),K_(x),&IN_ONE);}
}// GAIA


#endif
