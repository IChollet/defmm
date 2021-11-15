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

#ifndef GAIA_MATRICES_HPP
#define GAIA_MATRICES_HPP

#include <ostream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>

namespace gaia{

// Format Column-major pour les alias avec blas
template<typename T>
class Mat {
private:
  int N;
  int S;
  T* data;
public:
  Mat(){data = nullptr; N = 0; S = 0;}
  Mat(int _N, int _S) : N(_N) , S(_S) {data = new T[N*S];}
  Mat(int _N, int _S, const T &v) : N(_N) , S(_S) { data = new T[N*S];for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] = v;}}}
  Mat(int _N, int _S, const Mat &v)  : N(_N) , S(_S) { data = new T[N*S];for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] = v[i + j*N];}}}
  ~Mat(){if(data){std::free(data);} N = 0; S = 0;}
  void allocate(int _N, int _S){
    N = _N;
    S = _S;
    data = new T[N*S];
  }
  void allocate(int _N, int _S, T val){
    N = _N;
    S = _S;
    data = new T[N*S];
    for(int i = 0; i < N*S; i++){
      data[i] = val;
    }
  }
  void free(){
    if(data){std::free(data);}
    N = 0; S = 0;
  }
  const Mat &operator=(const T v) {for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] = v;}}return *this;}
  const Mat &operator+=(const T v) {for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] += v;}}return *this;}
  const Mat &operator-=(const T v) {for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] -= v;}}return *this;}
  const Mat &operator*=(const T v) {for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] *= v;}}return *this;}
  const Mat &operator/=(const T v) {for (int i=0; i<N; i++){for(int j=0;j<S;j++){ data[i + j*N] /= v;}}return *this;}
  const Mat &operator=(const Mat & v) {N = v.N; S = v.S; if(data){delete data;} data = new T[N*S];
    for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] = v[i + j*N];}}return *this;}
  const Mat &operator+=(const Mat & v) {for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] += v[i + j*N];}}return *this;}
  const Mat &operator-=(const Mat & v) {for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] -= v[i + j*N];}}return *this;}
  const Mat &operator*=(const Mat & v) {for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] *= v[i + j*N];}}return *this;}
  const Mat &operator/=(const Mat & v) {for(int i=0; i<N; i++){for(int j=0; j<S; j++){ data[i + j*N] /= v[i + j*N];}}return *this;}
  Mat operator+(const T v) const {return Mat(*this) += v;}
  Mat operator-(const T v) const {return Mat(*this) -= v;}
  Mat operator*(const T v) const {return Mat(*this) *= v;}
  Mat operator/(const T v) const {return Mat(*this) /= v;}
  Mat operator+(const Mat & v) const {return Mat(*this) += v;}
  Mat operator-(const Mat & v) const {return Mat(*this) -= v;}
  Mat operator*(const Mat & v) const {return Mat(*this) *= v;}
  Mat operator/(const Mat & v) const {return Mat(*this) /= v;}
  const T &operator()(int i, int j) const {return data[i + j*N];}
  T &operator()(int i, int j) {return data[i + j*N];}
  operator       T* ()       {return data;}
  operator const T* () const {return data;}
  friend T* K_(Mat& v){return v.data;}
  friend T* K_(Mat* v){return v->data;}
  friend std::ostream &operator<<(std::ostream & s, const Mat & v) {for (int i=0; i<v.N; i++){for(int j=0; j<v.S; j++){ s << v[i + j*v.N] << ' ';} s << std::endl;}return s;}
  friend inline int NbRow(const Mat & v){return v.N;}
  friend inline int NbCol(const Mat & v){return v.S;}
  friend inline int NbRow(const Mat * v){return v->N;}
  friend inline int NbCol(const Mat * v){return v->S;}
};
  
} // gaia

#endif
