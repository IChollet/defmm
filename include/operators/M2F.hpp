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

#ifndef M2F_MAIN_HEADER_HPP
#define M2F_MAIN_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{
  
  template<int DIM>
  inline void lfM2F(cell<DIM>& c){
    circulant::geom2fourier<DIM>(c.M, c.F);
  }

  template<int DIM>
  inline void hfM2F(cell<DIM>& c){
    for(int i = 0; i < c.ndir; i++){
      circulant::geom2fourier<DIM>(c.M+global::aLL*i,
				   c.F+global::aPP*i);
    }
  }
  
} // DEFMM


#endif
