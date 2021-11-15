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

#ifndef DUFMM_PLANEWAVE_HEADER_HPP
#define DUFMM_PLANEWAVE_HEADER_HPP

#include "../global.hpp"
#include "../thetis/thetis.hpp"
#include "../gaia/gaia.hpp"
#include "../structs/tree.hpp"

namespace defmm{

  namespace planewave{

    namespace internal{
      Vecc** Signatures;
      Vecc** SquaredSignatures;
    } // DEFMM::PLANEWAVE::INTERNAL

    template<int DIM>
    void computeSignatures(Tree<DIM>& T, Tree<DIM>& S){
      int nHFStages = 0;
      for(int l = 2; l < std::max(Nlvl(T),Nlvl(S)); l++){
        if(global::isHFStage[l]){nHFStages++;}
      }
      internal::Signatures        = new Vecc*[nHFStages];
      internal::SquaredSignatures = new Vecc*[nHFStages];
      for(int l = 2; l < std::max(Nlvl(T),Nlvl(S)); l++){
        if(global::precomputationAtStage[l] && global::isHFStage[l]){
	  internal::Signatures[global::HFStage[l]] = new Vecc[global::nDirPerStage[l]];
	  internal::SquaredSignatures[global::HFStage[l]] = new Vecc[global::nDirPerStage[l]];
	  flt locrad = Radius(T)/flt(1 << l);
	  for(int i = 0; i < global::nDirPerStage[l]; i++){
	    internal::Signatures[global::HFStage[l]][i].allocate(global::LL);
	    internal::SquaredSignatures[global::HFStage[l]][i].allocate(global::LL);
	    vecf<DIM> d = global::directions.tables[global::HFStage[l]][i].direction;
	    for(int k = 0; k < global::LL; k++){
	      veci<DIM> I = thetis::keys::i2I<DIM>(k,global::L);
	      vecf<DIM> x = thetis::uniform::point<DIM>(I,global::L,.5*locrad);
	      internal::Signatures[global::HFStage[l]][i][k] = exp(-cplxi*global::kappa*doprod(d,x));
	      internal::SquaredSignatures[global::HFStage[l]][i][k] = internal::Signatures[global::HFStage[l]][i][k] * internal::Signatures[global::HFStage[l]][i][k];
	    }
	  }
	}
      }
    }

    void applyNeg(int hfstg, int dir, cplx* in, cplx* out){
      cplx* e = K_(internal::Signatures[global::HFStage[hfstg]][dir]);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	out[i] = e[i] * in[i];
      }
    }

    void applyNeg2(int hfstg, int dir, cplx* in, cplx* out){
      cplx* e = K_(internal::SquaredSignatures[global::HFStage[hfstg]][dir]);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	out[i] = e[i] * in[i];
      }
    }

    void applyPos(int hfstg, int dir, cplx* in, cplx* out){
      cplx* e = K_(internal::Signatures[global::HFStage[hfstg]][dir]);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	out[i] = std::conj(e[i]) * in[i];
      }
    }

    void applyPos2(int hfstg, int dir, cplx* in, cplx* out){
      cplx* e = K_(internal::SquaredSignatures[global::HFStage[hfstg]][dir]);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	out[i] = std::conj(e[i]) * in[i];
      }
    }
    
    template<int DIM>
    inline void signature(cplx* in, cplx* out, vecf<DIM>& d, const flt& rad, const flt& sign){
      for(int i = 0; i < global::LL; i++){
	veci<DIM> I = thetis::keys::i2I<DIM>(i,global::L);
	vecf<DIM> x = thetis::uniform::point<DIM>(I,global::L,rad);
        out[i] = exp(sign*cplxi*global::kappa*doprod(d,x))*in[i];
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    template<int DIM>
    void computeDiffSignatures(Tree<DIM>& T, Tree<DIM>& S){
      int nHFStages = 0;
      for(int l = 2; l < std::max(Nlvl(T),Nlvl(S)); l++){
	if(global::precomputationAtStage[l] && global::isHFStage[l]){nHFStages++;}
      }
      internal::Signatures = new Vecc*[nHFStages];
      for(int l = 2; l < std::max(Nlvl(T),Nlvl(S)); l++){
        if(global::precomputationAtStage[l] && global::isHFStage[l] && !(global::HFStage[l]) ){
	  internal::Signatures[global::HFStage[l]] = new Vecc[global::nDirPerStage[l]];
	  flt locrad = Radius(T)/flt(1 << l);
	  for(int i = 0; i < global::nDirPerStage[l]; i++){
	    internal::Signatures[global::HFStage[l]][i].allocate(global::LL);
	    vecf<DIM> d = global::directions.tables[global::HFStage[l]][i].direction;
	    for(int k = 0; k < global::LL; k++){
	      veci<DIM> I = thetis::keys::i2I<DIM>(k,global::L);
	      vecf<DIM> x = thetis::uniform::point<DIM>(I,global::L,.5*locrad);
	      internal::Signatures[global::HFStage[l]][i][k] = exp(-cplxi*global::kappa*doprod(d,x));
	    }
	  }
	}
        if(global::precomputationAtStage[l] && global::isHFStage[l] && global::HFStage[l] ){
	  internal::Signatures[global::HFStage[l]] = new Vecc[global::nDirPerStage[l]];
	  flt locrad = Radius(T)/flt(1 << l);
	  for(int i = 0; i < global::nDirPerStage[l]; i++){
	    internal::Signatures[global::HFStage[l]][i].allocate(global::LL);
	    vecf<DIM> d  = global::directions.tables[global::HFStage[l]  ][i                         ].direction;
	    vecf<DIM> dt = global::directions.tables[global::HFStage[l-1]][int(floor(i/(1<<(DIM-1))))].direction;
	    vecf<DIM> d_ = d-dt;
	    for(int k = 0; k < global::LL; k++){
	      veci<DIM> I = thetis::keys::i2I<DIM>(k,global::L);
	      vecf<DIM> x = thetis::uniform::point<DIM>(I,global::L,.5*locrad);
	      internal::Signatures[global::HFStage[l]][i][k] = exp(-cplxi*global::kappa*doprod(d_,x));
	    }
	  }
	}
      }
    }
    
  } // PLANEWAVE
  
} // DEFMM

#endif
