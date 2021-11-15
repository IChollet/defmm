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

#ifndef THETIS_KEYS_HEADER_HPP
#define THETIS_KEYS_HEADER_HPP

#include "thetis_types.hpp"

namespace thetis{

  namespace keys{  
    
    template<int DIM>
    inline int I2i(const veci<DIM>& I, int L){
      int i = 0;
      for(int ii = 0; ii < DIM; ii++){
        i *= L;
        i += I[ii];
      }
      return i;
    }
    template<int DIM>
    inline veci<DIM> i2I(int i, int L){
      veci<DIM> I = 0;
      for(int ii=0; ii<DIM; ii++){
	I[DIM-ii-1] = (i%L);
	i -= I[DIM-ii-1];
	i /= L;
      }
      return I;
    }
    template<int DIM, int L>
    inline int I2i(const veci<DIM>& I){
      int i = 0;
      for(int ii = 0; ii < DIM; ii++){
	i *= L;
	i += I[ii];
      }
      return i;
    }
    template<int DIM, int L>
    inline veci<DIM> i2I(int i){
      veci<DIM> I = 0;
      for(int ii=0; ii<DIM; ii++){
	I[DIM-ii-1] = (i%L);
	i -= I[DIM-ii-1];
	i /= L;
      }
      return I;
    }
    template<int DIM>
    veci<DIM> NegPart(const veci<DIM>& In){
      veci<DIM> Out;
      for(int d = 0; d < DIM; d++){
	Out[d] = ( ( In[d] > -1 )  ?  0  :  -In[d] );
      }
      return Out;
    }
    template<int DIM>
    veci<DIM> PosPart(const veci<DIM>& In){
      veci<DIM> Out;
      for(int d = 0; d < DIM; d++){
	Out[d] = ( ( In[d] > -1 )  ?  In[d]  :  0 );
      }
      return Out;
    }
  } // keys
  
} // thetis


#endif
