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

#ifndef MAP_ON_TREES_HPP
#define MAP_ON_TREES_HPP

#include "../structs/tree.hpp"
#include "../calculus/planewave.hpp"
#include "../global.hpp"
#include <omp.h>

namespace defmm{

  template<int DIM>
  void fastComputeSourceLeafMatrices(Tree<DIM>& S){
    cell<DIM>* tab_cells = Root(S);
    for(int i = 0; i < nCells(S); i++){
      cell<DIM>* ci = tab_cells+i;
      if(!ci->son.size()){
	// Get Zmat
	Vecc** Z = new Vecc*[DIM];
	for(int k = 0; k < DIM; k++){
	  Z[k] = new Vecc[global::L];
	  for(int l = 0; l < global::L; l++){
	    Z[k][l].allocate(ci->nprt);
	  }
	}
        for(int k = 0; k < DIM; k++){
	  for(int l = 0; l < global::L; l++){
	    for(int j = 0; j < ci->nprt; j++){
	      Z[k][l][j] = thetis::uniform::S1D(l,ci->prt[j].pos[k]-ci->ctr[k],global::L,ci->rad);
	    }
	  }
	}
	// Assemble the matrix
        ci->Zmat.allocate(global::LL,ci->nprt);
	for(int l = 0; l < global::LL; l++){
	  veci<DIM> I = thetis::keys::i2I<DIM>(l,global::L);
	  for(int j = 0; j < ci->nprt; j++){
	    ci->Zmat(l,j) = 1.;
	    for(int k = 0; k < DIM; k++){
	      ci->Zmat(l,j) *= Z[k][I[k]][j];
	    }
	  }
	}
        if(global::isHFStage[ci->lvl] && global::precomputationAtStage[ci->lvl]){
	  for(int iD = 0; iD < global::nDirPerStage[ci->lvl]; iD++){
	    if(ci->effDir[iD]){
	      cplx* bufs  = new cplx[ci->nprt];
	      vecf<DIM> d = global::directions.tables[global::HFStage[ci->lvl]][iD].direction;
	      for(int j = 0; j < ci->nprt; j++){
		vecf<DIM> diff = ci->prt[j].pos-ci->ctr;
		bufs[j] = exp(-cplxi*global::kappa*doprod(d,diff));
	      }
	      ci->signatures.push_back(bufs);
	    }
	  }
	}
      }
    }
  }

  template<int DIM>
  void fastComputeTargetLeafMatrices(Tree<DIM>& T){
    cell<DIM>* tab_cells = Root(T);
    for(int i = 0; i < nCells(T); i++){
      cell<DIM>* ci = tab_cells+i;
      if(!ci->son.size()){
	// Allocate local table
	Vecc** Z = new Vecc*[DIM];
	for(int k = 0; k < DIM; k++){
	  Z[k] = new Vecc[global::L];
	  for(int l = 0; l < global::L; l++){
	    Z[k][l].allocate(ci->nprt);
	  }
	}
	// Map particles to the tensorized grid
	for(int k = 0; k < DIM; k++){
	  for(int l = 0; l < global::L; l++){
	    for(int i = 0; i < ci->nprt; i++){
	      Z[k][l][i] = thetis::uniform::S1D(l,ci->prt[i].pos[k] - ci->ctr[k],global::L,ci->rad);
	    }
	  }
	}
	ci->Zmat.allocate(ci->nprt,global::LL);
	for(int l = 0; l < global::LL; l++){
	  veci<DIM> I = thetis::keys::i2I<DIM>(l,global::L);
	  cplx* mL2P = K_(ci->Zmat)+l*ci->nprt;
#pragma omp simd
	  for(int j = 0; j < ci->nprt; j++){
	    mL2P[j] = 1.;
	    for(int k = 0; k < DIM; k++){
	      mL2P[j] *= Z[k][I[k]][j];
	    }
	  }
	}
        if(global::isHFStage[ci->lvl] && global::precomputationAtStage[ci->lvl]){
	  for(int iD = 0; iD < global::nDirPerStage[ci->lvl]; iD++){
	    if(ci->effDir[iD]){
	      cplx* bufs  = new cplx[ci->nprt];
	      vecf<DIM> d = global::directions.tables[global::HFStage[ci->lvl]][iD].direction;
	      for(int j = 0; j < ci->nprt; j++){
		vecf<DIM> diff = ci->prt[j].pos-ci->ctr;
		bufs[j] = exp(cplxi*global::kappa*doprod(d,diff));
	      }
	      ci->signatures.push_back(bufs);
	    }
	  }
	}
      }
    }
  }

  

} // DEFMM

#endif
