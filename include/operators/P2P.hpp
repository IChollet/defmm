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

#ifndef P2P_HEADER_HPP
#define P2P_HEADER_HPP

#include <omp.h>
#include "../structs/particles.hpp"
#include "../global.hpp"

namespace defmm{

  template<int DIM>
  void simdPartialDeinterleaveP2P(cell<DIM>& t, cell<DIM>& s){
    std::cout << f_red << "simdPartialDeinterleaveP2P has to be defined for DIM = " << DIM << f_def << std::endl;
    std::cout << f_red << "Exit(1)." << f_def << std::endl;
    exit(1);
  }

  template<>
  void simdPartialDeinterleaveP2P<3>(cell<3>& t, cell<3>& s){
    int N = t.nprt;
    int M = s.nprt;
    flt* pt;
    flt* ps;
    flt  pi4m1 = 1./pi4;
    flt x, y, z, R, K, r, co, si, pr, pi, qr, qi;
    cplx *pp = global::p+t.ind;
    cplx *qq = global::q+s.ind;
#pragma omp simd
    for(int i = 0; i < N; i++){
      pt = K_(t.prt[i].pos);
      pr = 0.;
      pi = 0.;
      for(int j = 0; j < M; j++){
	ps     = K_(s.prt[j].pos);
	qr     = qq[j].real();   // Deinterleave data
	qi     = qq[j].imag();
        x      = pt[0] - ps[0];
        y      = pt[1] - ps[1];
        z      = pt[2] - ps[2];
        R      = x*x+y*y+z*z;
	K      = 1./sqrt(R);
        if(R < 1.e-16){K = 0.;}
	r      = K * R;
	r     *= global::kappa;
	K     *= pi4m1;
	co     = cos(r);
	si     = sin(r);
	co    *= K;
	si    *= K;
	pr    += co*qr - si*qi;
	pi    += co*qi + si*qr;
      }
      pp[i] += cplx(pr,pi);
    }
  }
  
  template<>
  void simdPartialDeinterleaveP2P<2>(cell<2>& t, cell<2>& s){
    int N = t.nprt;
    int M = s.nprt;
    flt* pt;
    flt* ps;
    flt  pi4m1 = 1./pi4;
    flt x, y, R, K, r, co, si, pr, pi, qr, qi;
    cplx *pp = global::p+t.ind;
    cplx *qq = global::q+s.ind;
#pragma omp simd
    for(int i = 0; i < N; i++){
      pt = K_(t.prt[i].pos);
      pr = 0.;
      pi = 0.;
      for(int j = 0; j < M; j++){
	ps     = K_(s.prt[j].pos);
	qr     = qq[j].real();   // Deinterleave data
	qi     = qq[j].imag();
        x      = pt[0] - ps[0];
        y      = pt[1] - ps[1];
        R      = x*x+y*y;
	K      = 1./sqrt(R);     // Fast approx inverse sqrt
        if(R < 1.e-16){K = 0.;}
	r      = K * R;
	r     *= global::kappa;
	K     *= pi4m1;
	co     = cos(r);
	si     = sin(r);
	co    *= K;
	si    *= K;
	pr    += co*qr - si*qi;
	pi    += co*qi + si*qr;
      }
      pp[i] += cplx(pr,pi);
    }
  }
  template<>
  void simdPartialDeinterleaveP2P<1>(cell<1>& t, cell<1>& s){
    int N = t.nprt;
    int M = s.nprt;
    flt* pt;
    flt* ps;
    flt  pi4m1 = 1./pi4;
    flt x, R, K, r, co, si, pr, pi, qr, qi;
    cplx *pp = global::p+t.ind;
    cplx *qq = global::q+s.ind;
#pragma omp simd
    for(int i = 0; i < N; i++){
      pt = K_(t.prt[i].pos);
      pr = 0.;
      pi = 0.;
      for(int j = 0; j < M; j++){
	ps     = K_(s.prt[j].pos);
	qr     = qq[j].real();   // Deinterleave data
	qi     = qq[j].imag();
        x      = pt[0] - ps[0];
        R      = x*x;
	K      = 1./sqrt(R);     // Fast approx inverse sqrt
        if(R < 1.e-16){K = 0.;}
	r      = K * R;
	r     *= global::kappa;
	K     *= pi4m1;
	co     = cos(r);
	si     = sin(r);
	co    *= K;
	si    *= K;
	pr    += co*qr - si*qi;
	pi    += co*qi + si*qr;
      }
      pp[i] += cplx(pr,pi);
    }
  }
  
} // DEFMM

#endif
