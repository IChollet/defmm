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

#ifndef P2M_HEADER_HPP
#define P2M_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../calculus/planewave.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{
  
  // Tensorized version
  template<int DIM>
  inline void lfP2M(cell<DIM>& c){
    // Allocate local table
    Vecc** Z = new Vecc*[DIM];
    for(int k = 0; k < DIM; k++){
      Z[k] = new Vecc[global::L];
      for(int l = 0; l < global::L; l++){
	Z[k][l].allocate(c.nprt);
      }
    }
    // Map particles to the tensorized grid
    for(int k = 0; k < DIM; k++){
      for(int l = 0; l < global::L; l++){
#pragma omp simd
	for(int j = 0; j < c.nprt; j++){
	  Z[k][l][j] = thetis::uniform::S1D(l,c.prt[j].pos[k]-c.ctr[k],global::L,c.rad);
	}
      }
    }
    // Evaluate P2M using the tensorized structure
    for(int i = 0; i < global::LL; i++){
      veci<DIM> I = thetis::keys::i2I<DIM>(i,global::L);
      for(int j = 0; j < c.nprt; j++){
	cplx tmp = (global::q+c.ind)[j];
	for(int k = 0; k < DIM; k++){
	  tmp *= Z[k][I[k]][j];
	}
	c.M[i] += tmp;
      }
    }
  }

  template<int DIM>
  void hfP2M(cell<DIM>& c){
    // Allocate local table
    Vecc** Z = new Vecc*[DIM];
    for(int k = 0; k < DIM; k++){
      Z[k] = new Vecc[global::L];
      for(int l = 0; l < global::L; l++){
	Z[k][l].allocate(c.nprt);
      }
    }
    for(int iD = 0; iD < global::nDirPerStage[c.lvl]; iD++){
      if(c.effDir[iD]){
        vecf<DIM> d = global::directions.tables[global::HFStage[c.lvl]][iD].direction;
	// Map particles to the tensorized grid
	for(int k = 0; k < DIM; k++){
	  for(int l = 0; l < global::L; l++){
	    for(int j = 0; j < c.nprt; j++){
	      flt diffk = c.prt[j].pos[k]-c.ctr[k];
	      Z[k][l][j] = thetis::uniform::S1D(l,diffk,global::L,c.rad)
		* exp(-cplxi*global::kappa*d[k]*diffk);
	    }
	  }
	}
        // Evaluate P2M using the tensorized structure
        cplx* res = (c.M+(c.effDir[iD]-1)*global::aLL);
	for(int i = 0; i < global::LL; i++){
	  veci<DIM> I = thetis::keys::i2I<DIM>(i,global::L);
	  for(int j = 0; j < c.nprt; j++){
	    cplx tmp = (global::q+c.ind)[j];
	    for(int k = 0; k < DIM; k++){
	      tmp *= Z[k][I[k]][j];
	    }
	    res[i] += tmp;
	  }
	  res[i] *= exp(cplxi*global::kappa*doprod(d,
	  					   thetis::uniform::point<DIM>(I,
	  								       global::L,
	  								       c.rad)));
	}
      }
    }
  }

  /////////////////////////////
  template<int DIM>
  void phfP2M(cell<DIM>& c){
    Matc bufq(c.nprt,c.ndir);
    Matc bufp(global::LL,c.ndir);
    for(int k = 0; k < c.ndir; k++){
      cplx *bq = K_(bufq) + k*c.nprt;
#pragma omp simd
      for(int j = 0; j < c.nprt; j++){
        bq[j] = c.signatures[k][j]*(global::q+c.ind)[j];
      }
    }
    gaia::gemm(cplx(1.),K_(c.Zmat),K_(bufq),cplx(0.),K_(bufp),global::LL,c.nprt,c.ndir);
    for(int k = 0; k < c.ndir; k++){
      planewave::applyPos2(c.lvl,
			   c.dirs[k],
			   K_(bufp)+k*global::LL,
			   c.M+k*global::aLL);
    }
  }

  template<int DIM>
  void plfP2M(cell<DIM>& c){
    gaia::gemv(cplx(1.),K_(c.Zmat),global::q+c.ind,cplx(0.),c.M,global::LL,c.nprt);
  }
 
} // DEFMM


#endif
