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

#ifndef M2F_PASS_MAIN_HEADER_HPP
#define M2F_PASS_MAIN_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../operators/M2F.hpp"
#include "../structs/initdevs.hpp"
#include <omp.h>

namespace defmm{

  template<int DIM>
  inline void M2FPass(cell<DIM>* Cells, int nCells){
    for(int c = 0; c < nCells; c++){
      if(global::precomputationAtStage[Cells[c].lvl]){
	if(global::isHFStage[Cells[c].lvl]){
	  hfM2F(Cells[c]);
        }else{
	  lfM2F(Cells[c]);
        }
      }
    }
  }

  template<int DIM>
  void M2FPass(Tree<DIM>& T){
#ifdef DEFMM_VERBOSE
    flt T0,T1;
    T0 = gaia::chronos::time();
#endif
    M2FPass(Root(T),nCells(T));
#ifdef DEFMM_VERBOSE
    T1 = gaia::chronos::time();
    std::cout << t_bold << f_yellow << "M2F pass: " << T1-T0  << t_def << f_def << std::endl;
#endif
  }

  ////////////////////////////////////////////////////////////////

  namespace internal{
    fftw_plan P_frwrd;
  } // INTERNAL
  
  template<int DIM> void InitM2F(){
    int* Parray = new int[DIM];
    for(int k = 0; k < DIM; k++){Parray[k] = global::P;}
    internal::P_frwrd = fftw_plan_many_dft(DIM,
					 Parray,
					 internal::nSrcDevs,
					 reinterpret_cast<fftw_complex*>(global::Fs),
					 NULL,
					 1,
					 global::aPP,
					 reinterpret_cast<fftw_complex*>(global::Fs),
					 NULL,
					 1,
					 global::aPP,
					 FFTW_FORWARD,
					 FFTW_ESTIMATE);
  }
  
  template<int DIM>
  void M2FPassBATCH(){
    flt t0,t1,t2;
    t0 = gaia::chronos::time();
    for(int p = 0; p < internal::nSrcDevs; p++){
      cplx* intmp = global::Fs + p*global::aPP;
      cplx* in    = global::Ms + p*global::aLL;
      for(int i = 0; i < global::aPP; i++){
        intmp[i] = cplx0;
      }
      for(int i = 0; i < global::LL; i++){
        intmp[circulant::internal::itilda[i]] = in[i];
      }
    }
    t1 = gaia::chronos::time();
    fftw_execute_dft(internal::P_frwrd,
		     reinterpret_cast<fftw_complex*>(global::Fs),
		     reinterpret_cast<fftw_complex*>(global::Fs));
    t2 = gaia::chronos::time();
    std::cout << t_bold << f_red << t1-t0 << "\t" << t2-t1 << t_def << f_def << std::endl;
  }
  
} // DEFMM

#endif
