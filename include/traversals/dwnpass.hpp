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

#ifndef DOWNWARD_PASS_HEADER_HPP
#define DOWNWARD_PASS_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../operators/L2P.hpp"
#include "../operators/L2L.hpp"

namespace defmm{

  template<int DIM>
  void DownwardPass(int c, cell<DIM>* Cells){
    if (isLeaf(Cells[c])){
      if(global::precomputationAtStage[(Cells+c)->lvl]){
	if(global::isHFStage[(Cells+c)->lvl]){
	  phfL2P(Cells[c]);
        }else{
	  plfL2P(Cells[c]);
	}
      }
    }else{
      for(unsigned int i = 0; i < Cells[c].son.size(); i++){
        if(global::precomputationAtStage[(Cells+c)->lvl]){
	  if(global::isHFStage[(Cells+c)->lvl]){
	    if(global::HFStage[(Cells+c)->lvl]){
	      hfL2Lblock(Cells[Cells[c].son[i]], Cells[c]);
	    }else{
	      hf2lfL2L(Cells[Cells[c].son[i]], Cells[c]);
	    }
	  }else{
	    lfL2L(Cells[Cells[c].son[i]], Cells[c]);
	  }
	}
	DownwardPass(Cells[c].son[i],Cells);
      }
    }
  }

  template<int DIM>
  void DownwardPass(Tree<DIM>& T){
#ifdef DEFMM_VERBOSE
    flt t0,t1;
    t0 = gaia::chronos::time();
#endif
    DownwardPass(0,Root(T));
#ifdef DEFMM_VERBOSE
    t1 = gaia::chronos::time();
    std::cout << f_yellow << t_bold << "Dwn pass: " << t1-t0 << t_def << f_def << std::endl;
#endif
  }

  //////////////////////////////////////////////////////////////////////////////////////
  ////// BLANK DWN PASS

  template<int DIM>
  void BlankDownwardPass(int c, cell<DIM>* Cells){
    // Quit if c is a leaf
    if(!Cells[c].son.size()){return;}
    // Propagate the direction information
    int ns = global::directions.nsons[global::HFStage[Cells[c].lvl+1]];
    if(global::isHFStage[Cells[c].lvl] && global::precomputationAtStage[Cells[c].lvl]){
      for(unsigned int i = 0; i < Cells[c].son.size(); i++){
	if(global::isHFStage[Cells[c].lvl+1]){
	  for(int dt = 0; dt < global::nDirPerStage[Cells[c].lvl+1]; dt++){
	    for(int k = 0; k < ns; k++){
	      if(Cells[c].effDir[ns*dt+k]){
		Cells[Cells[c].son[i]].effDir[dt] = 1;
	      }
	    }
	  }
	  BlankDownwardPass<DIM>(Cells[c].son[i],Cells);
	}
      }
    }else{
      for(unsigned int i = 0; i < Cells[c].son.size(); i++){
	BlankDownwardPass<DIM>(Cells[c].son[i],Cells);
      }
    }
  }

  template<int DIM>
  void BlankDownwardPass(Tree<DIM>& T){
    BlankDownwardPass<DIM>(0,Root(T));
  }
  
  
} // DEFMM

#endif
