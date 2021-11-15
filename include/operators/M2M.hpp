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

#ifndef M2M_HEADER_HPP
#define M2M_HEADER_HPP

#include "../calculus/interpolation.hpp"
#include "../symmetries/hypercube.hpp"

namespace defmm{
  
  template<int DIM>
  inline void lfM2M(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = s.ctr - t.ctr;
    Matc bufq(global::LL,1);
    fastlagrange::loadBuf (s.M,K_(bufq));
    fastlagrange::realAnterpolation<DIM>(bufq,v);
    fastlagrange::writeBuf(t.M,K_(bufq));
  }
  
  
  template<int DIM>
  void hfM2M(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = s.ctr - t.ctr;
    int tlvl    = t.lvl;
    int slvl    = s.lvl;
    int thflvl  = global::HFStage[tlvl];
    int shflvl  = global::HFStage[slvl];
    int ns      = global::directions.nsons[shflvl];
    Vecc bufq(global::LL);
    for(int ds = 0; ds < global::nDirPerStage[s.lvl]; ds++){
      for(int dt = 0; dt < ns; dt++){
	int dirind = ds*ns+dt;
	if(t.effDir[dirind]){
	  vecf<DIM> d = global::directions.tables[thflvl][dirind].direction;
	  planewave::applyNeg(tlvl,
			      dirind,
			      s.M+(s.effDir[ds]-1)*global::aLL,
			      K_(bufq));
	  bufq *= exp(-cplxi*global::kappa*doprod(d,v));
	  fastlagrange::anterpolation<DIM>(bufq,v);
	  planewave::applyPos2(tlvl,dirind,K_(bufq),K_(bufq));
	  cplx* res = (t.M+(t.effDir[dirind]-1)*global::aLL);
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    res[i] += bufq[i];
	  }
	}
      }
    }
  }

  template<int DIM>
  void hfM2Mblock(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = s.ctr - t.ctr;
    int tlvl    = t.lvl;
    int slvl    = s.lvl;
    int thflvl  = global::HFStage[tlvl];
    int shflvl  = global::HFStage[slvl];
    int ns      = global::directions.nsons[shflvl];
    Matc BUFQ(global::LL,t.ndir);
    for(int k = 0; k < t.ndir; k++){
      planewave::applyNeg(tlvl,
			  t.dirs[k],
			  s.M+(s.effDir[t.dirs[k]/ns]-1)*global::aLL,
			  K_(BUFQ)+k*global::LL);
    }
    fastlagrange::realAnterpolation<DIM>(BUFQ,v);
    for(int k = 0; k < t.ndir; k++){
      cplx* u   = K_(BUFQ)+global::LL*k;
      planewave::applyPos2(tlvl,
			   t.dirs[k],
			   u,
			   u);
      vecf<DIM> d = global::directions.tables[thflvl][t.dirs[k]].direction;
      cplx edv = exp(-cplxi*global::kappa*doprod(d,v));
      cplx* res = (t.M+(t.effDir[t.dirs[k]]-1)*global::aLL);
#pragma omp simd
      for(int i = 0; i < global::LL; i++){
	res[i] += u[i]*edv;
      }
    }
  }

  template<int DIM>
  void lf2hfM2M(cell<DIM>& t, cell<DIM>& s){
    vecf<DIM> v = s.ctr - t.ctr;
    Vecc bufq(global::LL);
    for(int dt = 0; dt < global::nDirPerStage[t.lvl]; dt++){
      if(t.effDir[dt]){
	vecf<DIM> d = global::directions.tables[global::HFStage[t.lvl]][dt].direction;
	planewave::applyNeg(t.lvl,dt,s.M,K_(bufq));
	bufq *= exp(-cplxi*global::kappa*doprod(d,v));
	fastlagrange::anterpolation<DIM>(bufq,v);
        planewave::applyPos2(t.lvl,dt,bufq,bufq);
	cplx* res = (t.M+(t.effDir[dt]-1)*global::aLL);
	for(int i = 0; i < global::LL; i++){
	  res[i] += bufq[i];
	}
      }
    }
  }

} // DEFMM

#endif
