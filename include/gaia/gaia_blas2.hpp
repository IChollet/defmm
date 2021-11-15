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

#ifndef GAIA_BLAS_2_HPP
#define GAIA_BLAS_2_HPP

#include <complex>
#include "gaia_vectors.hpp"
#include "gaia_matrices.hpp"
#include "gaia_blas_consts.hpp"

namespace gaia{
  extern "C"{
    void sgemv_(const char*, int*, int*, BLAS_S*,BLAS_S*, int*, BLAS_S*, int*, BLAS_S*, BLAS_S*, int*);
    void dgemv_(const char*, int*, int*, BLAS_D*,BLAS_D*, int*, BLAS_D*, int*, BLAS_D*, BLAS_D*, int*);
    void cgemv_(const char*, int*, int*, BLAS_C*,BLAS_C*, int*, BLAS_C*, int*, BLAS_C*, BLAS_C*, int*);
    void zgemv_(const char*, int*, int*, BLAS_Z*,BLAS_Z*, int*, BLAS_Z*, int*, BLAS_Z*, BLAS_Z*, int*);
  }
  inline void gemv(BLAS_S u, Mat<BLAS_S>& A, Vec<BLAS_S>& x, BLAS_S v, Vec<BLAS_S>& y){int n = NbRow(A); int s = NbCol(A); sgemv_(charN, &n, &s, &u, K_(A), &n, K_(x), &IN_ONE, &v, K_(y), &IN_ONE);}
  inline void gemv(BLAS_D u, Mat<BLAS_D>& A, Vec<BLAS_D>& x, BLAS_D v, Vec<BLAS_D>& y){int n = NbRow(A); int s = NbCol(A); dgemv_(charN, &n, &s, &u, K_(A), &n, K_(x), &IN_ONE, &v, K_(y), &IN_ONE);}
  inline void gemv(BLAS_C u, Mat<BLAS_C>& A, Vec<BLAS_C>& x, BLAS_C v, Vec<BLAS_C>& y){int n = NbRow(A); int s = NbCol(A); cgemv_(charN, &n, &s, &u, K_(A), &n, K_(x), &IN_ONE, &v, K_(y), &IN_ONE);}
  inline void gemv(BLAS_Z u, Mat<BLAS_Z>& A, Vec<BLAS_Z>& x, BLAS_Z v, Vec<BLAS_Z>& y){int n = NbRow(A); int s = NbCol(A); zgemv_(charN, &n, &s, &u, K_(A), &n, K_(x), &IN_ONE, &v, K_(y), &IN_ONE);}

  inline void gemv(Mat<BLAS_S>& A, Vec<BLAS_S>& x, Vec<BLAS_S>& y){int n = NbRow(A); int s = NbCol(A); sgemv_(charN, &n, &s, &S_ONE, K_(A), &n, K_(x), &IN_ONE, &S_ZERO, K_(y), &IN_ONE);}
  inline void gemv(Mat<BLAS_D>& A, Vec<BLAS_D>& x, Vec<BLAS_D>& y){int n = NbRow(A); int s = NbCol(A); dgemv_(charN, &n, &s, &D_ONE, K_(A), &n, K_(x), &IN_ONE, &D_ZERO, K_(y), &IN_ONE);}
  inline void gemv(Mat<BLAS_C>& A, Vec<BLAS_C>& x, Vec<BLAS_C>& y){int n = NbRow(A); int s = NbCol(A); cgemv_(charN, &n, &s, &C_ONE, K_(A), &n, K_(x), &IN_ONE, &C_ZERO, K_(y), &IN_ONE);}
  inline void gemv(Mat<BLAS_Z>& A, Vec<BLAS_Z>& x, Vec<BLAS_Z>& y){int n = NbRow(A); int s = NbCol(A); zgemv_(charN, &n, &s, &Z_ONE, K_(A), &n, K_(x), &IN_ONE, &Z_ZERO, K_(y), &IN_ONE);}

  inline void gemv(BLAS_S* A, BLAS_S* x, BLAS_S* y, int L){int n = L; int s = L; sgemv_(charN, &n, &s, &S_ONE, A, &n, x, &IN_ONE, &S_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_D* A, BLAS_D* x, BLAS_D* y, int L){int n = L; int s = L; dgemv_(charN, &n, &s, &D_ONE, A, &n, x, &IN_ONE, &D_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_C* A, BLAS_C* x, BLAS_C* y, int L){int n = L; int s = L; cgemv_(charN, &n, &s, &C_ONE, A, &n, x, &IN_ONE, &C_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_Z* A, BLAS_Z* x, BLAS_Z* y, int L){int n = L; int s = L; zgemv_(charN, &n, &s, &Z_ONE, A, &n, x, &IN_ONE, &Z_ZERO, y, &IN_ONE);}

  inline void gemv(BLAS_S* A, BLAS_S* x, BLAS_S* y, int L, int K){int n = L; int s = K; sgemv_(charN, &n, &s, &S_ONE, A, &n, x, &IN_ONE, &S_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_D* A, BLAS_D* x, BLAS_D* y, int L, int K){int n = L; int s = K; dgemv_(charN, &n, &s, &D_ONE, A, &n, x, &IN_ONE, &D_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_C* A, BLAS_C* x, BLAS_C* y, int L, int K){int n = L; int s = K; cgemv_(charN, &n, &s, &C_ONE, A, &n, x, &IN_ONE, &C_ZERO, y, &IN_ONE);}
  inline void gemv(BLAS_Z* A, BLAS_Z* x, BLAS_Z* y, int L, int K){int n = L; int s = K; zgemv_(charN, &n, &s, &Z_ONE, A, &n, x, &IN_ONE, &Z_ZERO, y, &IN_ONE);}

  inline void gemv(BLAS_S u, BLAS_S* A, BLAS_S* x, BLAS_S v, BLAS_S* y, int L, int K){int n = L; int s = K; sgemv_(charN, &n, &s, &u, A, &n, x, &IN_ONE, &v, y, &IN_ONE);}
  inline void gemv(BLAS_D u, BLAS_D* A, BLAS_D* x, BLAS_D v, BLAS_D* y, int L, int K){int n = L; int s = K; dgemv_(charN, &n, &s, &u, A, &n, x, &IN_ONE, &v, y, &IN_ONE);}
  inline void gemv(BLAS_C u, BLAS_C* A, BLAS_C* x, BLAS_C v, BLAS_C* y, int L, int K){int n = L; int s = K; cgemv_(charN, &n, &s, &u, A, &n, x, &IN_ONE, &v, y, &IN_ONE);}
  inline void gemv(BLAS_Z u, BLAS_Z* A, BLAS_Z* x, BLAS_Z v, BLAS_Z* y, int L, int K){int n = L; int s = K; zgemv_(charN, &n, &s, &u, A, &n, x, &IN_ONE, &v, y, &IN_ONE);}

  namespace algo{
    // Produit restrictif :
    // sb : Indice du premier élément de B ciblé
    // sc : Indice du premier élément de C à considérer
    // Cette fonction effectue le produit matrice/vecteur d'une partie de B par A dans une partie de C
    inline void gemv(Mat<BLAS_S>& A, Vec<BLAS_S>& B, int sb, Vec<BLAS_S>& C, int sc){int n = NbRow(A); int s = NbCol(A);sgemv_(charN, &n, &s, &S_ONE, K_(A), &n, K_(B)+sb, &IN_ONE, &S_ZERO, K_(C)+sc, &IN_ONE);}
    inline void gemv(Mat<BLAS_D>& A, Vec<BLAS_D>& B, int sb, Vec<BLAS_D>& C, int sc){int n = NbRow(A); int s = NbCol(A);dgemv_(charN, &n, &s, &D_ONE, K_(A), &n, K_(B)+sb, &IN_ONE, &D_ZERO, K_(C)+sc, &IN_ONE);}
    inline void gemv(Mat<BLAS_C>& A, Vec<BLAS_C>& B, int sb, Vec<BLAS_C>& C, int sc){int n = NbRow(A); int s = NbCol(A);cgemv_(charN, &n, &s, &C_ONE, K_(A), &n, K_(B)+sb, &IN_ONE, &C_ZERO, K_(C)+sc, &IN_ONE);}
    inline void gemv(Mat<BLAS_Z>& A, Vec<BLAS_Z>& B, int sb, Vec<BLAS_Z>& C, int sc){int n = NbRow(A); int s = NbCol(A);zgemv_(charN, &n, &s, &Z_ONE, K_(A), &n, K_(B)+sb, &IN_ONE, &Z_ZERO, K_(C)+sc, &IN_ONE);}
  } // ALGO
} // GAIA

#endif
