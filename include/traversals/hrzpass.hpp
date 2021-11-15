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

#ifndef HORIZONTAL_PASS_HEADER_HPP
#define HORIZONTAL_PASS_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../operators/M2L.hpp"
#include "../operators/P2P.hpp"

namespace defmm{

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// ADMISSIBILITY CONDITIONS
  
  template<int DIM>
  inline bool MAClist(cell<DIM>& X, cell<DIM>& Y){
    veci<DIM> I = X.idx - Y.idx;
    int test  = 0;
    int ncmp  = DIM;
    for(int k = 0; k < DIM; k++){
      ncmp -= (I[k] ? 0 : 1);
      test += std::abs(I[k]);
    }
    return (test > ncmp);
  }

  template<int DIM>
  inline bool MACgeomhf(cell<DIM>& X, cell<DIM>& Y){
    flt tmp = (1.-1.e-8)*(norm2(X.ctr - Y.ctr)-2.*Y.rad);
    flt tt = 2.*Y.rad;
    return (tmp >= std::max(tt*tt*global::kappa, tt));
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// HORIZONTAL PASS
  
  template<int DIM>
  void HorizontalPass(int i, int j, cell<DIM>* X, cell<DIM>* Y){
    if(
       global::precomputationAtStage[X[i].lvl] &&
       X[i].lvl == Y[j].lvl
       ){
      if((X+i)->lvl && global::isHFStage[(X+i)->lvl]){
        if(MACgeomhf(X[i], Y[j])){
	  hfM2L<DIM>(X[i],Y[j]);
	  return;
	}
      }else{
     	if(MAClist(X[i], Y[j])){
	  lfM2L<DIM>(X[i],Y[j]);
	  return;
	}
      }
    }
    if(
        isLeaf(X[i]) ||
        isLeaf(Y[j])
       ){
      simdPartialDeinterleaveP2P(X[i],Y[j]);
      return;
    }
    for(unsigned int k1 = 0; k1 < X[i].son.size(); k1++){
      for(unsigned int k2 = 0; k2 < Y[j].son.size(); k2++){
	HorizontalPass(X[i].son[k1],Y[j].son[k2],X,Y);
      }
    }
  }
  
  template<int DIM>
  void HorizontalPass(Tree<DIM>& T, Tree<DIM>& S){
#ifdef DEFMM_VERBOSE
    flt t0,t1;
    t0 = gaia::chronos::time();
#endif
    HorizontalPass(0,0,Root(T),Root(S));
#ifdef DEFMM_VERBOSE
    t1 = gaia::chronos::time();
    std::cout << t_bold << f_yellow << "Hrz pass: " << t1-t0 << f_def << t_def << std::endl;
#endif
  }


  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// HF-M2L BLANK HORIZONTAL PASS

  template<int DIM>
  void BlankHorizontalPass(int i, int j, cell<DIM>* X, cell<DIM>* Y){
    // Quit the pass in the low frequency regime
    if(
       (!global::isHFStage[X[i].lvl] && !global::isHFStage[Y[j].lvl]) ||
       (!global::isHFStage[X[i].lvl] && !Y[j].son.size()) ||
       (!X[i].son.size() && !global::isHFStage[Y[j].lvl])
       ){
      return;
    }
    // Find interaction if in HF regime
    if(
       (X[i].lvl == Y[j].lvl) &&
       MACgeomhf(X[i], Y[j])  &&
       global::precomputationAtStage[X[i].lvl]
       ){
      blankM2Lhf(X[i], Y[j]);
      return;
    }
    // Recursion
    if(X[i].son.size() && Y[j].son.size()){
      for(unsigned int k1 = 0; k1 < X[i].son.size(); k1++){
        for(unsigned int k2 = 0; k2 < Y[j].son.size(); k2++){
	  BlankHorizontalPass(X[i].son[k1],Y[j].son[k2],X,Y);
	}
      }
    }
  }
  
  template<int DIM>
  void BlankHorizontalPass(Tree<DIM>& T, Tree<DIM>& S){
    BlankHorizontalPass(0,0,Root(T),Root(S));
  }
  
} // DEFMM

#endif
