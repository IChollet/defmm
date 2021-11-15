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

#ifndef L2L_HEADER_HPP
#define L2L_HEADER_HPP

#include "../calculus/interpolation.hpp"
#include "../symmetries/hypercube.hpp"

namespace defmm{

  template<int DIM>
  inline void lfL2L(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = t.ctr - s.ctr;
    Vecc bufq(global::LL);
    fastlagrange::loadBuf (s.L,bufq);
    fastlagrange::interpolation<DIM>(bufq,v);
    fastlagrange::writeBuf(t.L,bufq);
  }

  template<int DIM>
  void hfL2L(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = t.ctr - s.ctr;
    int tlvl    = t.lvl;
    int slvl    = s.lvl;
    int thflvl  = global::HFStage[tlvl];
    int shflvl  = global::HFStage[slvl];
    int ns      = global::directions.nsons[thflvl];
    Vecc bufq(global::LL);
    for(int dt = 0; dt < global::nDirPerStage[tlvl]; dt++){
      for(int ds = 0; ds < ns; ds++){
	int dirind = dt*ns+ds;
	if(s.effDir[dirind]){
	  vecf<DIM> d = global::directions.tables[shflvl][dirind].direction;
	  planewave::applyNeg2(slvl,
			       dirind,
			       s.L+(s.effDir[dirind]-1)*global::aLL,
			       K_(bufq));
	  fastlagrange::interpolation<DIM>(bufq,v);
	  planewave::applyPos(slvl,
			      dirind,
			      K_(bufq),
			      K_(bufq));
	  cplx ev = exp(cplxi*global::kappa*doprod(v,d));
	  cplx* res = (t.L+(t.effDir[dt]-1)*global::aLL);
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    res[i] += bufq[i]*ev;
	  }
	}
      }
    }
  }

  template<int DIM>
  void hfL2Lblock(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = t.ctr - s.ctr;
    int tlvl    = t.lvl;
    int slvl    = s.lvl;
    int thflvl  = global::HFStage[tlvl];
    int shflvl  = global::HFStage[slvl];
    int ns      = global::directions.nsons[thflvl];
    Matc BUFQ(global::LL,s.ndir);
    for(int k = 0; k < s.ndir; k++){
      planewave::applyNeg2(slvl,
			   s.dirs[k],
			   s.L+(s.effDir[s.dirs[k]]-1)*global::aLL,
			   K_(BUFQ)+global::LL*k);
    }
    fastlagrange::realInterpolation<DIM>(BUFQ,v);
    //fastlagrange::interpolation<DIM>(BUFQ,v);
    for(int k = 0; k < s.ndir; k++){
      cplx* u   = K_(BUFQ)+global::LL*k;
      planewave::applyPos(slvl,
			  s.dirs[k],
			  u,
			  u);
      vecf<DIM> d = global::directions.tables[shflvl][s.dirs[k]].direction;
      cplx  edv = exp(cplxi*global::kappa*doprod(v,d));
      cplx* res = (t.L+(t.effDir[s.dirs[k]/ns]-1)*global::aLL);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	res[i] += u[i]*edv;
      }
    }
  }

  template<int DIM>
  void hf2lfL2L(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = t.ctr - s.ctr;
    Vecc bufq(global::LL);
    for(int ds = 0; ds < global::nDirPerStage[s.lvl]; ds++){
      if(s.effDir[ds]){
	vecf<DIM> d = global::directions.tables[global::HFStage[s.lvl]][ds].direction;
	planewave::applyNeg2(s.lvl,ds,s.L+(s.effDir[ds]-1)*global::aLL,K_(bufq));
        fastlagrange::interpolation<DIM>(bufq,v);
	planewave::applyPos(s.lvl,ds,K_(bufq),K_(bufq));
	cplx ev = exp(cplxi*global::kappa*doprod(v,d));
#pragma omp simd
	for(int i = 0; i < global::LL; i++){
	  t.L[i] += bufq[i]*ev;
	}
      }
    }
  }
  
} // DEFMM

#endif
