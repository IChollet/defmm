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

#ifndef UPWARD_PASS_HEADER_HPP
#define UPWARD_PASS_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../operators/P2M.hpp"
#include "../operators/M2M.hpp"

namespace defmm{

  template<int DIM>
  void UpwardPass(int c, cell<DIM>* Cells){
    if(isLeaf(Cells[c])){
      if(global::precomputationAtStage[(Cells+c)->lvl]){
	if(global::isHFStage[(Cells+c)->lvl]){
	  phfP2M(Cells[c]);
        }else{
	  plfP2M(Cells[c]);
        }
      }
    }else{
      for(unsigned int i = 0; i < Cells[c].son.size(); i++){
        UpwardPass(Cells[c].son[i],Cells);
      }
      for(unsigned int i = 0; i < Cells[c].son.size(); i++){
        if(global::precomputationAtStage[(Cells+c)->lvl]){
	  if(global::isHFStage[(Cells+c)->lvl]){
	    if(global::HFStage[(Cells+c)->lvl]){
	      hfM2Mblock(Cells[c], Cells[Cells[c].son[i]]);
	    }else{
	      lf2hfM2M(Cells[c], Cells[Cells[c].son[i]]);
	    }
	  }else{
	    lfM2M(Cells[c], Cells[Cells[c].son[i]]);
	  }
	}
      }
    }
  }

  template<int DIM>
  void UpwardPass(Tree<DIM>& T){
#ifdef DEFMM_VERBOSE
    flt T0,T1;
    T0 = gaia::chronos::time();
#endif
    UpwardPass(0,Root(T));
#ifdef DEFMM_VERBOSE
    T1 = gaia::chronos::time();
    std::cout << f_yellow << t_bold << "UPW pass: " << T1-T0 << f_def << t_def << std::endl;
#endif
  }
  
} // DEFMM

#endif
