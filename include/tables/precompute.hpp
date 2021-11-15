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

#ifndef DUFMM_PRECOMPUTE_HEADER_HPP
#define DUFMM_PRECOMPUTE_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../symmetries/permutations.hpp"
#include "frequency.hpp"
#include "../symmetries/hypercube.hpp"

namespace defmm{

  template<int DIM, thetis::kernel_op K>
  void prcmputeM2LTable(Tree<DIM>& T, Tree<DIM>& S){
    int nlvl = std::max(Nlvl(T),Nlvl(S));
    for(int l = 2; l < nlvl; l++){
      if(global::precomputationAtStage[l] && !global::isHFStage[l]){
	flt locrad = Radius(T)/flt(1 << l);
	int cpt = 0;
	for(int p = 0; p < mypow(global::nCellPerDirMax[l],DIM); p++){
          bool b = false;
	  veci<DIM> iV = thetis::keys::i2I<DIM>(p,global::nCellPerDirMax[l]);
	  for(int k = 0; k < DIM; k++){b = (b || (iV[k] >= global::nCellPerDirMin[l]));}
	  for(int k = 0; k < DIM-1; k++){b = (b && (iV[k] >= iV[k+1]));}
	  if(b){
	    cpt++;
	    vecf<DIM> fV;
	    for(int k = 0; k < DIM; k++){
	      fV[k] = 2.*flt(iV[k])*locrad;
	    }
#ifdef M2L_LF_NO_INDIRECTIONS
	    Vecc tmpVecc(global::PP); tmpVecc = cplx0;
	    int *permut;
	    circulant::fillCircCoeffs<DIM,K>(tmpVecc, locrad, fV);
	    for(int per = 0; per < (1 << (DIM*(DIM-1)/2+DIM)); per++){
	      permut = K_(permutations::table[per]);
#pragma omp simd
	      for(int j = 0; j < global::PP; j++){
		K_(global::lf_M2L_tables[l][p][per])[j] = K_(tmpVecc)[permutations::table[per][j]];
	      }
	    }
#else
	    circulant::fillCircCoeffs<DIM,K>(global::M2LTable[l][p], locrad, fV);
#endif
	  }
	}
#ifdef DEFMM_VERBOSE
	std::cout << "nb precomputed M2L at level " << f_cyan << l << f_def << " : " << f_cyan << cpt << f_def << std::endl;
#endif
      }
    }
  }

} // DEFMM

#endif 
