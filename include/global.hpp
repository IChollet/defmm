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
//  (see ./../LICENSE.txt)
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with defmm.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================

#ifndef GLOBAL_VARIABLES_MAIN_HEADER_HPP
#define GLOBAL_VARIABLES_MAIN_HEADER_HPP

#include "./tables/newdirections.hpp"

//#define INITIALIZE_DEV_BLOCK
#define M2L_CSTE_INDIRECTIONS
//#define M2L_LF_NO_INDIRECTIONS
#define DEFMM_VERBOSE


namespace defmm{

  inline int mypow(int a, int b){
    int tmp = 1;
    for(int i = 1; i <= b; i++){tmp *= a;}
    return tmp;
  }

  namespace global{

    // Problem parameters
    flt kappa;    // Wavenumber
    int L;        // 1D interpolation order
    int P;        // 2*L-1 -> 1D padding
    int LL;       // pow(L,DIM) -> number of pts
                  // in the geometric tensorized
                  // grid.
    int PP;       // pow(P,DIM) -> number of pts
                  // in the tensorized fourier grid.
    int aLL;      // aLL > LL and nextelem / alignof(cplx) \in \bb{N}.
    int aPP;      // Same for PP. (needed for padding)

    // Inner FMM vectors
    cplx* p;
    cplx* q;

    // Expansion containers
    cplx* Ms;
    cplx* Fs;
    cplx* Ls;
    cplx* Js;

    // Global tables
    DirectionTables<3> directions;     // Tree and tables of directions
#ifdef M2L_CSTE_INDIRECTIONS
    Vecc**             M2LTable;       // Where to store precomputed M2Ls
#else
    Vecc***            M2LTable;       // Where to store precomputed M2Ls
#endif
#ifdef M2L_LF_NO_INDIRECTIONS
    Vecc***            lf_M2L_tables;
#endif
    int***             dirTable;       // Each M2L corresponds to a given direction
                                       // (in the fundamental domain)
    int*               nCellPerDirMax; // Maximal distance (along a basis
                                       // axis) between two admissible cells
                                       // for each stages.
    int*               nCellPerDirMin; // Minimal distance
    int*               isHFStage;      // True if a stage is HF
    int*               HFStage;        // Effective HF stage
    bool*              precomputationAtStage; // is there any precomp at stage ...
    int*               nDirPerStage;   // Nb direction per stage
    bool               noHF;           // true if there is no HF regime
    int                nHFStage;       // Number of HF stages

    // Parameters setting
    template<int DIM>
    void setParameters(int _L, flt _kappa){
      L        = _L;
      kappa    = _kappa;
      LL       = mypow(L,DIM);
      P        = 2*L-1;
      PP       = mypow(P,DIM);     
      aLL      = LL + ((LL%alignof(cplx)) ? alignof(cplx)-(LL%alignof(cplx)) : 0);
      aPP      = PP + ((PP%alignof(cplx)) ? alignof(cplx)-(PP%alignof(cplx)) : 0);
    }
    
    
  } // GLOBAL

} // DEFMM


#endif
