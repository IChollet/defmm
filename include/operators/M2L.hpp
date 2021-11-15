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

#ifndef M2L_HEADER_HPP
#define M2L_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../thetis/thetis.hpp"
#include "../calculus/circulant.hpp"
#include "../tables/precompute.hpp"
#include <omp.h>

namespace defmm{
  
  /*
    Allocate memory for M2L tables on fundamental domain.
  */
  template<int DIM, typename M2L_TYPE>
  void allocateM2LTable(Tree<DIM>& T, Tree<DIM>& S){
    int nlvl         = std::max(Nlvl(T), Nlvl(S));
#ifdef M2L_LF_NO_INDIRECTIONS
    global::lf_M2L_tables = new Vecc**[nlvl];
#endif
    global::M2LTable = new M2L_TYPE*[nlvl];
    global::dirTable = new int**[nlvl];
    global::nHFStage = 0;
    for(int l = 0; l < nlvl; l++){
      global::nHFStage += (global::isHFStage[l] ? 1 : 0);
    }
#ifdef M2L_LF_NO_INDIRECTIONS
    int nPer = (1 << (DIM+DIM*(DIM-1)/2));
    for(int l = std::max(global::nHFStage, 2); l < nlvl; l++){
      long int nPlh = mypow(global::nCellPerDirMax[l],DIM);
      global::lf_M2L_tables[l] = new Vecc*[nPlh];
      for(int p = 0; p < nPlh; p++){
	global::lf_M2L_tables[l][p] = new Vecc[nPer];
	for(int q = 0; q < nPer; q++){
	  global::lf_M2L_tables[l][p][q].allocate(global::PP,64);
	  global::lf_M2L_tables[l][p][q] = cplx0;
	}
      }
    }
#endif
    for(int l = 1; l < nlvl; l++){
      if(global::precomputationAtStage[l]){
	int nPl = mypow(4,DIM);
	if(global::isHFStage[l] || global::isHFStage[l-1]){
	  nPl = mypow(global::nCellPerDirMax[l],DIM);
	}
	global::M2LTable[l] = new M2L_TYPE[nPl];
	global::dirTable[l] = new int*[nPl];
	if(!global::isHFStage[l]){
	  for(int p = 0; p < nPl; p++){
	    global::M2LTable[l][p].allocate(global::PP,64);
	    global::M2LTable[l][p] = cplx0;
	  }
	}
      }
    }
  }
  
  /*
    Get effective M2L and return index of associated permutation.
  */
  template<int DIM, typename M2L_EFF_TYPE>
  int* chooseM2L(int lvl, veci<DIM>& iZ, M2L_EFF_TYPE* effM2L){
    int *permutation;
    veci<DIM> nI;
    getPermuteAndFundamentalMultiIndex(iZ,nI,permutation);
    int iz = thetis::keys::I2i<DIM>(nI,global::nCellPerDirMax[lvl]);
    *effM2L = K_(global::M2LTable[lvl][iz]);
    return permutation;
  }

  /*
    Get effective M2L and return index of associated permutation.
    Also gives the index of best cone (iobc).
  */
  template<int DIM, typename M2L_EFF_TYPE>
  int* chooseM2Lhf(int lvl, veci<DIM>& iZ, M2L_EFF_TYPE* effM2L, int* iobc){
    int *permutation;
    veci<DIM> nI;
    int hsh;
    getPermuteAndFundamentalMultiIndexAndHashcode(iZ,nI,permutation,hsh);
    int iz = thetis::keys::I2i<DIM>(nI,global::nCellPerDirMax[lvl]);
    *effM2L = K_(global::M2LTable[lvl][iz]);
    *iobc   = global::dirTable[lvl][iz][hsh];
    return permutation;
  }
  
  template<int DIM>
  void hfM2L(cell<DIM>& t, cell<DIM>& s){
    veci<DIM> iZ = t.idx - s.idx;
    int iobc;
    cplx *u;
    int  *p = chooseM2Lhf(t.lvl, iZ, &u,&iobc);
    cplx *y = s.F+(s.effDir[iobc]-1)*global::aPP;
    cplx *x = t.J+(t.effDir[iobc]-1)*global::aPP;
#pragma omp simd
    for(int i = 0; i < global::PP; i++){
      x[i] +=  u[p[i]]*y[i];
    }
  }
  
  template<int DIM>
  void lfM2L(cell<DIM>& t, cell<DIM>& s){
#ifdef M2L_LF_NO_INDIRECTIONS
    veci<DIM> iZ = t.idx - s.idx;
    int hashcode0 = hashcodeZov2Z(iZ);
    int hashcode1 = hashcodeSigma(iZ);
    int hashcode2 = hashcode0 + (hashcode1<<DIM);
    int* per = permutations::sigma[hashcode1];
    veci<DIM> R;
    for(int k = 0; k < DIM; k++){R[k] = std::abs(iZ[per[k]]);}
    int  iz = thetis::keys::I2i<DIM>(R,global::nCellPerDirMax[t.lvl]);
    cplx *u = K_(global::lf_M2L_tables[t.lvl][iz][hashcode2]);
    cplx *y = s.F;
    cplx *x = t.J;
#pragma omp simd
    for(int i = 0; i < global::PP; i++){
      x[i] += u[i]*y[i];
    }
#else
    veci<DIM> iZ = t.idx - s.idx;
    cplx *u;
    int  *p = chooseM2L(t.lvl, iZ, &u);
    cplx *y = s.F;
    cplx *x = t.J;
#pragma omp simd
    for(int i = 0; i < global::PP; i++){
      x[i] += u[p[i]]*y[i];
    }
#endif
  }
  
  namespace hrzpassinternal{
    int cpt;
  } // HRZPASSINTERNAL
  
  template<int DIM>
  void blankM2Lhf(cell<DIM>& t, cell<DIM>& s){
    // Permute the translation vector
    veci<DIM> iZ = t.idx - s.idx;
    int hsh = hashcodeZov2Z(iZ) + (hashcodeSigma(iZ)<<DIM);
    veci<DIM> nI;
    int *permutation;
    getPermuteAndFundamentalMultiIndex(iZ,nI,permutation);
    int iz = thetis::keys::I2i<DIM>(nI,global::nCellPerDirMax[t.lvl]);
    vecf<DIM> fV;
    flt locrad = t.rad;
    for(int k = 0; k < DIM; k++){
      fV[k] = 2.*flt(nI[k])*locrad;
    }
    // Precompute the corresponding M2L
    if(!global::M2LTable[t.lvl][iz]){
      global::M2LTable[t.lvl][iz].allocate(global::PP);
      global::M2LTable[t.lvl][iz] = cplx0;
      circulant::fillCircCoeffs<DIM,thetis::kernels::Helmholtz>(global::M2LTable[t.lvl][iz], locrad, fV);
      global::dirTable[t.lvl][iz] = new int[permutations::nPerm];
      for(int ll = 0; ll < permutations::nPerm; ll++){
	global::dirTable[t.lvl][iz][ll] = -1;
      }
      hrzpassinternal::cpt++;
    }
    // Add direction to the source and the target cells
    if(global::dirTable[t.lvl][iz][hsh] == -1){
      int iobc;
      vecf<DIM> fZ = t.ctr - s.ctr;
      fZ /= norm2(fZ);
      iobc = bestConeIndex<DIM>(global::directions,
				global::HFStage[t.lvl],
				fZ);
      global::dirTable[t.lvl][iz][hsh] = iobc;
    }
    t.effDir[global::dirTable[t.lvl][iz][hsh]] = 1;
    s.effDir[global::dirTable[t.lvl][iz][hsh]] = 1;
  }
  
} // DEFMM


#endif
