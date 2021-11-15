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

#ifndef F2L_MAIN_HEADER_HPP
#define F2L_MAIN_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{
  
  template<int DIM>
  inline void lfF2L(cell<DIM>& c){
    circulant::fourier2geom<DIM>(c.J, c.L);
  }
  
  template<int DIM>
  inline void hfF2L(cell<DIM>& c){
    for(int i = 0; i < c.ndir; i++){
      circulant::fourier2geom<DIM>(c.J+global::aPP*i,
				   c.L+global::aLL*i);
    }
  }

} // DEFMM


#endif
