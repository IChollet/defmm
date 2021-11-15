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

#ifndef CIRCULANT_HEADER_HPP
#define CIRCULANT_HEADER_HPP

#include "../global.hpp"
#include "../thetis/thetis.hpp"
#include "../thetis/titan.hpp"
#include <omp.h>

namespace defmm{

  namespace circulant{

    namespace internal{
      int* itilda;
      flt* tabZ;
      fftw_plan p_forward;
      fftw_plan p_backward;

    } // INTERNAL
    
    template<int DIM>
    void precomputeCirculant(){
      int Parray[DIM];
      for(int k = 0; k < DIM; k++){Parray[k] = global::P;}
      cplx* intest = NULL;
      cplx* outest = NULL;
      if(posix_memalign((void**)&intest, 64, sizeof(cplx)) != 0){
	std::cerr << "Allocation failed (in 'precomputeCirculant' in circulant.hpp)" << std::endl;
	exit(1);
      }
      if(posix_memalign((void**)&outest, 64, sizeof(cplx)) != 0){
	std::cerr << "Allocation failed (in 'precomputeCirculant' in circulant.hpp)" << std::endl;
	exit(1);
      }
      internal::p_forward  = fftw_plan_dft(DIM,
					   Parray,
					   reinterpret_cast<fftw_complex*>(intest),
					   reinterpret_cast<fftw_complex*>(outest),
					   FFTW_FORWARD,
					   FFTW_ESTIMATE);
      internal::p_backward = fftw_plan_dft(DIM,
					   Parray,
					   reinterpret_cast<fftw_complex*>(intest),
					   reinterpret_cast<fftw_complex*>(outest),
					   FFTW_BACKWARD,
					   FFTW_ESTIMATE);
      
      internal::itilda = new int[global::PP];
      for(int i = 0; i < global::LL; i++){
	internal::itilda[i] = thetis::keys::I2i<DIM>(thetis::keys::i2I<DIM>(i, global::L), global::P);
      }
      int size = global::PP*DIM;
      internal::tabZ = new flt[size];
      int indInTabZ = 0;
      for(int i = 0; i < global::PP; i++){
	veci<DIM> E = thetis::keys::i2I<DIM>(i,global::P);
        for(int k = 0; k < DIM; k++){
	  E[k] = ( (E[k] < global::L) ? E[k] : -(2*global::L-1 - E[k]) );
	}
	veci<DIM> I = thetis::keys::PosPart<DIM>(E);
	veci<DIM> J = thetis::keys::NegPart<DIM>(E);
	vecf<DIM> X = thetis::uniform::point<DIM>(I,global::L,1.);
	vecf<DIM> Y = thetis::uniform::point<DIM>(J,global::L,1.);
        for(int k = 0; k < DIM; k++){
	  internal::tabZ[indInTabZ++] = X[k]-Y[k];
	}
      }
    }
    
    template<int DIM>
    inline void geom2fourier(cplx* in, cplx* out){
      Vecc intmp(global::PP);
      intmp = cplx(0.);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
        intmp[internal::itilda[i]] = in[i];
      }
      fftw_execute_dft(internal::p_forward,
		       reinterpret_cast<fftw_complex*>(K_(intmp)),
		       reinterpret_cast<fftw_complex*>(out));
    }
    
    template<int DIM>
    inline void fourier2geom(cplx* in, cplx* out){
      Vecc intmp(global::PP);
      fftw_execute_dft(internal::p_backward,
		       reinterpret_cast<fftw_complex*>(in),
		       reinterpret_cast<fftw_complex*>(K_(intmp)));
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
        out[i] += intmp[internal::itilda[i]]/flt(global::PP);
      }
    }

    template<int DIM, thetis::kernel_op K>
    inline void fillCircCoeffs(Vecc& out, const flt& rad, vecf<DIM>& v){
      Vecc in(global::PP);
      int indInTabZ = 0;
      for(int i = 0; i < global::PP; i++){
	vecf<DIM> Z;
	for(int k = 0; k < DIM; k++){Z[k] = internal::tabZ[indInTabZ++]*rad;}
        flt R = norm2(v+Z);
	in[i] = K(R,global::kappa);
      }
      titan::fft::forward<DIM>(K_(in),K_(out),global::P);
    }

    template<>
    void fillCircCoeffs<3,thetis::kernels::Helmholtz>(Vecc& out, const flt& rad, vecf<3>& v){
      Vecc in(global::PP);
      int indInZ = 0;
      for(int i = 0; i < global::PP; i++){
	vecf<3> Z;
	for(int k = 0; k < 3; k++){Z[k] = internal::tabZ[indInZ++]*rad;}
        flt R = norm2(v+Z);
	in[i] = exp(cplxi*R*global::kappa)/(pi4*R);
      }
      titan::fft::forward<3>(K_(in),K_(out),global::P);
    }
    
  } // CIRCULANT
  
} // DEFMM

#endif
