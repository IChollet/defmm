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

#ifndef PARTICLES_AND_CELLS_HEADER_HPP
#define PARTICLES_AND_CELLS_HEADER_HPP

#include <iostream>
#include <complex>
#include <vector>
#include <stdint.h>
#include <fstream>
#include <string>
#include "../gaia/gaia.hpp"

namespace defmm{

  template <int DIM>
  struct particle{
    vecf<DIM> pos; // Position
    int       ind; // Correspondance in matrix / vector entries
  };
  template <int DIM> using particles = std::vector<particle<DIM> >;
  
  template <int DIM>
  struct cell{
    veci<DIM>          idx;          // Integer coordinates in the stage
    int                lvl;          // Level of the cell
    int                nprt;         // Number of particles
    vecf<DIM>          ctr;          // Center
    flt                rad;          // Radius
    std::vector<int>   son;          // Indices of sons in the cell array
    particle<DIM>*     prt;          // Pointer on first particle
    cplx*              M;            // Multipole expansions
    cplx*              F;            // Fourier multipole expansions
    cplx*              J;            // Fourier local expansions
    cplx*              L;            // Local expansions
    Veci               effDir;       // Effective directions
    int                ndir;         // Nombre de directions effectives
    int*               dirs;         // Array of ndir effective direction indices
    std::vector<cplx*> signatures;   // Signatures at hf-leaf level
    Matc               Zmat;         // Interpolation matrix at leaves
    int                ind; // Index in the FMM numerotation of first particle/first term in vectors
    
    std::vector<cell*> P2PSrcs;
    std::vector<cplx*> P2PMats;

  };
  template <int DIM> using cells = std::vector<cell<DIM> >;

  template<int DIM>
  inline bool isLeaf(cell<DIM>& c){return (c.son.size()==0);}

  /*
    This function code is based on the function 'getBounds'
    in exafmm (https://github.com/exafmm/exafmm/blob/master/include/build_tree.h)
   */
  template <int DIM>
  void getBounds(particles<DIM>& tparts,
		 vecf<DIM>& tBoxCenter,
		 flt& tBoxRadius,
		 particles<DIM>& sparts,
		 vecf<DIM>& sBoxCenter,
		 flt& sBoxRadius){
    vecf<DIM> totmn = tparts[0].pos;
    vecf<DIM> totmx = tparts[0].pos;
    for(unsigned int i = 1; i < tparts.size(); i++){
      totmn = min(totmn,tparts[i].pos);
      totmx = max(totmx,tparts[i].pos);
    }
    for(unsigned int i = 0; i < sparts.size(); i++){
      totmn = min(totmn,sparts[i].pos);
      totmx = max(totmx,sparts[i].pos);
    }
    tBoxCenter = (totmn+totmx)/2;
    tBoxRadius = std::max(max(tBoxCenter-totmn), max(totmx-tBoxCenter));
    sBoxCenter = tBoxCenter;
    sBoxRadius = tBoxRadius;
  }
  
} // DEFMM

#endif
