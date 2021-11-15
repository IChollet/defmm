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

#ifndef L2P_HEADER_HPP
#define L2P_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{

  // Version tensorisée
  template<int DIM>
  void lfL2P(cell<DIM>& c){
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
	for(int i = 0; i < c.nprt; i++){
	  Z[k][l][i] = thetis::uniform::S1D(l,c.prt[i].pos[k] - c.ctr[k],global::L,c.rad);
	}
      }
    }
    // Evaluate L2P using the tensorized structure
    for(int j = 0; j < global::LL; j++){
      veci<DIM> J = thetis::keys::i2I<DIM>(j,global::L);
      for(int i = 0; i < c.nprt; i++){
	cplx tmp = c.L[j];
	for(int k = 0; k < DIM; k++){
	  tmp *= Z[k][J[k]][i];
	}
        (global::p+c.ind)[i] += tmp;
      }
    }
  }

  template<int DIM>
  void hfL2P(cell<DIM>& c){
    // Allocate local table
    Vecc** Z = new Vecc*[DIM];
    for(int k = 0; k < DIM; k++){
      Z[k] = new Vecc[global::L];
      for(int l = 0; l < global::L; l++){
	Z[k][l].allocate(c.nprt);
      }
    }
    for(int iD = 0; iD < global::nDirPerStage[c.lvl]; iD++){
      Vecc buff(global::LL);
      if(c.effDir[iD]){
	planewave::applyNeg(c.lvl,
			    iD,
			    K_(c.L)+(c.effDir[iD]-1)*global::aLL,
			    K_(buff));
	planewave::applyNeg(c.lvl,
			    iD,
			    K_(buff),
			    K_(buff));
       
	vecf<DIM> d = global::directions.tables[global::HFStage[c.lvl]][iD].direction;
	// Map particles to the tensorized grid
	// NOTE : CES ETAPES NE DEPENDENT PAS DE LA DIRECTION!!!
	for(int k = 0; k < DIM; k++){
	  for(int l = 0; l < global::L; l++){
	    for(int i = 0; i < c.nprt; i++){
	      Z[k][l][i] = thetis::uniform::S1D(l,c.prt[i].pos[k] - c.ctr[k],global::L,c.rad);
	    }
	  }
	}
	// Init a buffer
	Vecc bufp(c.nprt); bufp = cplx0;
	// Evaluate L2P using the tensorized structure
	for(int j = 0; j < global::LL; j++){
	  veci<DIM> J = thetis::keys::i2I<DIM>(j,global::L);
	  for(int i = 0; i < c.nprt; i++){
	    cplx tmp = buff[j];//c.L[j+iD*global::LL];
	    for(int k = 0; k < DIM; k++){
	      tmp *= Z[k][J[k]][i];
	    }
	    bufp[i] += tmp;
	  }
	}
	for(int i = 0; i < c.nprt; i++){
	  vecf<DIM> diff = c.prt[i].pos - c.ctr;
	  (global::p+c.ind)[i] += bufp[i] * exp( cplxi*global::kappa*doprod(d,diff));
	}
      }
    }
  }

  ///////////////////////////////////////////////////
  template<int DIM>
  void phfL2P(cell<DIM>& c){
    for(int iD = 0; iD < global::nDirPerStage[c.lvl]; iD++){
      Vecc buff(global::LL);
      if(c.effDir[iD]){
        planewave::applyNeg2(c.lvl,
			    iD,
			    c.L+(c.effDir[iD]-1)*global::aLL,
			    K_(buff));
	vecf<DIM> d = global::directions.tables[global::HFStage[c.lvl]][iD].direction;
	Vecc bufp(c.nprt);
	gemv(c.Zmat,buff,bufp);
	//gaia::gemv(1.,K_(c.Zmat),K_(buff),0.,K_(bufp),c.nprt,global::LL);
	///*
	for(int i = 0; i < c.nprt; i++){
	  vecf<DIM> diff = c.prt[i].pos - c.ctr;
	  (global::p+c.ind)[i] += bufp[i] * exp( cplxi*global::kappa*doprod(d,diff));
	}
	//*/
      }
    }
  }

  template<int DIM>
  void pphfL2P(cell<DIM>& c){
    Matc buff(global::LL,c.ndir);
    Matc bufp(c.nprt,c.ndir);
    for(int k = 0; k < c.ndir; k++){
      planewave::applyNeg2(c.lvl,
			   c.dirs[k],
			   c.L+(c.effDir[c.dirs[k]]-1)*global::aLL,
			   K_(buff)+k*global::LL);
    }
    gemm(c.Zmat,buff,bufp);
    for(int k = 0; k < c.ndir; k++){
      vecf<DIM> d = global::directions.tables[global::HFStage[c.lvl]][c.effDir[k]].direction;
      for(int i = 0; i < c.nprt; i++){
	vecf<DIM> diff = c.prt[i].pos - c.ctr;
	(global::p+c.ind)[i] += (K_(bufp)+k*global::LL)[i] * exp( cplxi*global::kappa*doprod(d,diff));
      }
    }
  }

  template<int DIM>
  void plfL2P(cell<DIM>& c){
    gaia::gemv(cplx(1.),
	       K_(c.Zmat),
	       c.L,
	       cplx(1.),
	       global::p+c.ind,
	       c.nprt,
	       global::LL);
  }

} // DEFMM


#endif
