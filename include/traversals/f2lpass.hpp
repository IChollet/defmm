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

#ifndef F2L_PASS_MAIN_HEADER_HPP
#define F2L_PASS_MAIN_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../operators/F2L.hpp"

namespace defmm{

  template<int DIM>
  void F2LPass(cell<DIM>* Cells, int nCells){
    for(int c = 0; c < nCells; c++){
      if(global::precomputationAtStage[Cells[c].lvl]){
        if(global::isHFStage[Cells[c].lvl]){
	  hfF2L<DIM>(Cells[c]);
        }else{
	  lfF2L<DIM>(Cells[c]);
        }
      }
    }
  }

  template<int DIM>
  void F2LPass(Tree<DIM>& T){
#ifdef DEFMM_VERBOSE
    flt t0,t1;
    t0 = gaia::chronos::time();
#endif
    F2LPass(Root(T),nCells(T));
#ifdef DEFMM_VERBOSE
    t1 = gaia::chronos::time();
    std::cout << f_yellow << t_bold << "F2L pass: " << t1-t0 << t_def << f_def << std::endl;
#endif
  }

} // DEFMM

#endif
