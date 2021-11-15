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

#ifndef GAIA_BLAS_3_HPP
#define GAIA_BLAS_3_HPP

#include <complex>
#include "gaia_vectors.hpp"
#include "gaia_matrices.hpp"
#include "gaia_blas_consts.hpp"

namespace gaia{
  extern "C"{
  void sgemm_(const char*, const char*, unsigned*, unsigned*, unsigned*, BLAS_S*, BLAS_S*, unsigned*, BLAS_S*, unsigned*, BLAS_S*, BLAS_S*, unsigned*);
  void dgemm_(const char*, const char*, unsigned*, unsigned*, unsigned*, BLAS_D*, BLAS_D*, unsigned*, BLAS_D*, unsigned*, BLAS_D*, BLAS_D*, unsigned*);
  void cgemm_(const char*, const char*, unsigned*, unsigned*, unsigned*, BLAS_C*, BLAS_C*, unsigned*, BLAS_C*, unsigned*, BLAS_C*, BLAS_C*, unsigned*);
  void zgemm_(const char*, const char*, unsigned*, unsigned*, unsigned*, BLAS_Z*, BLAS_Z*, unsigned*, BLAS_Z*, unsigned*, BLAS_Z*, BLAS_Z*, unsigned*);
  }
  inline void gemm(BLAS_S u, Mat<BLAS_S>& A, Mat<BLAS_S>& B, BLAS_S v, Mat<BLAS_S>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    sgemm_(charN,charN,&n,&s,&k,&u,K_(A),&n,K_(B),&k, &v,K_(C),&n);}
  inline void gemm(BLAS_D u, Mat<BLAS_D>& A, Mat<BLAS_D>& B, BLAS_D v, Mat<BLAS_D>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    dgemm_(charN,charN,&n,&s,&k,&u,K_(A),&n,K_(B),&k, &v,K_(C),&n);}
  inline void gemm(BLAS_C u, Mat<BLAS_C >& A, Mat<BLAS_C >& B, BLAS_C v, Mat<BLAS_C >& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    cgemm_(charN,charN,&n,&s,&k,&u,K_(A),&n,K_(B),&k, &v,K_(C),&n);}
  inline void gemm(BLAS_Z u, Mat<BLAS_Z >& A, Mat<BLAS_Z >& B, BLAS_Z v, Mat<BLAS_Z >& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    zgemm_(charN,charN,&n,&s,&k,&u,K_(A),&n,K_(B),&k, &v,K_(C),&n);}
  inline void gemm(Mat<BLAS_S>& A, Mat<BLAS_S>& B, Mat<BLAS_S>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    sgemm_(charN,charN,&n,&s,&k,&S_ONE,K_(A),&n,K_(B),&k, &S_ZERO,K_(C),&n);}
  inline void gemm(Mat<BLAS_D>& A, Mat<BLAS_D>& B, Mat<BLAS_D>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    dgemm_(charN,charN,&n,&s,&k,&D_ONE,K_(A),&n,K_(B),&k, &D_ZERO,K_(C),&n);}
  inline void gemm(Mat<BLAS_C>& A, Mat<BLAS_C>& B, Mat<BLAS_C>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    cgemm_(charN,charN,&n,&s,&k,&C_ONE,K_(A),&n,K_(B),&k, &C_ZERO,K_(C),&n);}
  inline void gemm(Mat<BLAS_Z>& A, Mat<BLAS_Z>& B, Mat<BLAS_Z>& C){
    unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = NbCol(B);
    zgemm_(charN,charN,&n,&s,&k,&Z_ONE,K_(A),&n,K_(B),&k, &Z_ZERO,K_(C),&n);}

  inline void gemm(BLAS_S u, BLAS_S* A, BLAS_S* B, BLAS_S v, BLAS_S* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    sgemm_(charN,charN,&n,&s,&k,&u,A,&n,B,&k, &v,C,&n);}
  inline void gemm(BLAS_D u, BLAS_D* A, BLAS_D* B, BLAS_D v, BLAS_D* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    dgemm_(charN,charN,&n,&s,&k,&u,A,&n,B,&k, &v,C,&n);}
  inline void gemm(BLAS_C u, BLAS_C* A, BLAS_C* B, BLAS_C v, BLAS_C* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    cgemm_(charN,charN,&n,&s,&k,&u,A,&n,B,&k, &v,C,&n);}
  inline void gemm(BLAS_Z u, BLAS_Z* A, BLAS_Z* B, BLAS_Z v, BLAS_Z* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    zgemm_(charN,charN,&n,&s,&k,&u,A,&n,B,&k, &v,C,&n);}
  inline void gemm(BLAS_S* A, BLAS_S* B, BLAS_S* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    sgemm_(charN,charN,&n,&s,&k,&S_ONE,A,&n,B,&k, &S_ZERO,C,&n);}
  inline void gemm(BLAS_D* A, BLAS_D* B, BLAS_D* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    dgemm_(charN,charN,&n,&s,&k,&D_ONE,A,&n,B,&k, &D_ZERO,C,&n);}
  inline void gemm(BLAS_C* A, BLAS_C* B, BLAS_C* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    cgemm_(charN,charN,&n,&s,&k,&C_ONE,A,&n,B,&k, &C_ZERO,C,&n);}
  inline void gemm(BLAS_Z* A, BLAS_Z* B, BLAS_Z* C, int _n, int _k, int _s){
    unsigned int n = _n;  unsigned int k = _k;  unsigned int s = _s;
    zgemm_(charN,charN,&n,&s,&k,&Z_ONE,A,&n,B,&k, &Z_ZERO,C,&n);}
  
  namespace algo{
    // gvmm suppose que B a un nombre d'entée mutiple du nombre de colones de A.
    // Il s'agit d'appliquer A aux Size(B)/n sous-vecteurs de B. 
    inline void gvmm(Mat<BLAS_Z >& A, Vec<BLAS_Z >& B, Vec<BLAS_Z >& C){
      unsigned int n = NbRow(A);  unsigned int k = NbCol(A);  unsigned int s = Size(B)/n;
      zgemm_(charN,charN,&n,&s,&k,&Z_ONE,K_(A),&n,K_(B),&k, &Z_ZERO,K_(C),&n);}
  }
  
} // GAIA

#endif
