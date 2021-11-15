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

#ifndef INITDEVS_HPP
#define INITDEVS_HPP

namespace defmm{
  
  template<int DIM>
  void InitializeSourceArrays(Tree<DIM>& S){
    for(int i = 0; i < nCells(S); i++){
      if(global::precomputationAtStage[(*Cells(S))[i].lvl]){
	if(global::isHFStage[(*Cells(S))[i].lvl]){
	  (*Cells(S))[i].effDir.allocate(global::nDirPerStage[(*Cells(S))[i].lvl]);
	  (*Cells(S))[i].effDir = 0;
        }
      }
    }
  }
  
  template<int DIM>
  void InitializeTargetArrays(Tree<DIM>& T){
    for(int i = 0; i < nCells(T); i++){
      if(global::precomputationAtStage[(*Cells(T))[i].lvl]){
	if(global::isHFStage[(*Cells(T))[i].lvl]){
	  (*Cells(T))[i].effDir.allocate(global::nDirPerStage[(*Cells(T))[i].lvl]);
	  (*Cells(T))[i].effDir = 0;
	}
      }
    }
  }

  
  namespace internal{
    int nSrcDevs;
    long int MsSize;
    long int FsSize;
  } // INTERNAL
  
  template<int DIM>
  void InitializeSourceDevsBlock(Tree<DIM>& S){
    // Count number of source dev
    internal::nSrcDevs = 0;
    for(int i = 0; i < nCells(S); i++){
      if(global::precomputationAtStage[(*Cells(S))[i].lvl]){
	if(global::isHFStage[(*Cells(S))[i].lvl]){
	  for(int k = 0; k < global::nDirPerStage[(*Cells(S))[i].lvl]; k++){
	    if((*Cells(S))[i].effDir[k]){
	      internal::nSrcDevs++;
	    }
	  }
        }else{
	  internal::nSrcDevs++;
	}
      }
    }
    // Initialize containers
    internal::MsSize = internal::nSrcDevs*global::aLL;
    internal::FsSize = internal::nSrcDevs*global::aPP;
    if(internal::MsSize < 0 || internal::FsSize < 0){std::cout << "InitializeSourceDevs (precompute.hpp) -> negative array size" << std::endl; exit(1);}
    global::Ms = new cplx[internal::MsSize];
    global::Fs = new cplx[internal::FsSize];
    // Associate pointers to cells and get direction indices
    int count = 0;
    for(int i = 0; i < nCells(S); i++){
      if(global::precomputationAtStage[(*Cells(S))[i].lvl]){
	if(global::isHFStage[(*Cells(S))[i].lvl]){
	  int cpt = 0;
	  for(int k = 0; k < global::nDirPerStage[(*Cells(S))[i].lvl]; k++){
	    if((*Cells(S))[i].effDir[k]){
	      (*Cells(S))[i].effDir[k] = ++cpt;
	    }
	  }
	  (*Cells(S))[i].ndir = cpt;
	  (*Cells(S))[i].dirs = new int[(*Cells(S))[i].ndir];
	  int loccount = 0;
	  for(int u = 0; u < global::nDirPerStage[(*Cells(S))[i].lvl]; u++){
	    if((*Cells(S))[i].effDir[u]){(*Cells(S))[i].dirs[loccount++] = u;}
	  }
	  (*Cells(S))[i].M    = global::Ms + count*global::aLL;
	  (*Cells(S))[i].F    = global::Fs + count*global::aPP;
	  count              += cpt;
        }else{
	  (*Cells(S))[i].M    = global::Ms + count*global::aLL;
	  (*Cells(S))[i].F    = global::Fs + count*global::aPP;
	  count++;
	}
      }
    }
    for(int i = 0; i < internal::MsSize; i++){
      global::Ms[i] = cplx0;
    }
    for(int i = 0; i < internal::FsSize; i++){
      global::Fs[i] = cplx0;
    }
  }
  
  
  namespace internal{
    int nTrgDevs;
    long int LsSize;
    long int JsSize;
  } // INTERNAL
  template<int DIM>
  void InitializeTargetDevsBlock(Tree<DIM>& T){
    // Count number of source dev
    internal::nTrgDevs = 0;
    for(int i = 0; i < nCells(T); i++){
      if(global::precomputationAtStage[(*Cells(T))[i].lvl]){
	if(global::isHFStage[(*Cells(T))[i].lvl]){
	  for(int k = 0; k < global::nDirPerStage[(*Cells(T))[i].lvl]; k++){
	    if((*Cells(T))[i].effDir[k]){
	      internal::nTrgDevs++;
	    }
	  }
        }else{
	  internal::nTrgDevs++;
	}
      }
    }
    // Initialize containers
    internal::LsSize = internal::nTrgDevs*global::aLL;
    internal::JsSize = internal::nTrgDevs*global::aPP;
    if(internal::LsSize < 0 || internal::JsSize < 0){std::cout << "InitializeTargetDevs (precompute.hpp) -> negative array size" << std::endl; exit(1);}
    global::Ls = new cplx[internal::LsSize];
    global::Js = new cplx[internal::JsSize];
    // Associate pointers to cells and get direction indices
    int count = 0;
    for(int i = 0; i < nCells(T); i++){
      if(global::precomputationAtStage[(*Cells(T))[i].lvl]){
	if(global::isHFStage[(*Cells(T))[i].lvl]){
	  int cpt = 0;
	  for(int k = 0; k < global::nDirPerStage[(*Cells(T))[i].lvl]; k++){
	    if((*Cells(T))[i].effDir[k]){
	      (*Cells(T))[i].effDir[k] = ++cpt;
	    }
	  }
	  (*Cells(T))[i].ndir = cpt;
	  (*Cells(T))[i].dirs = new int[(*Cells(T))[i].ndir];
	  int loccount = 0;
	  for(int u = 0; u < global::nDirPerStage[(*Cells(T))[i].lvl]; u++){
	    if((*Cells(T))[i].effDir[u]){(*Cells(T))[i].dirs[loccount++] = u;}
	  }
	  (*Cells(T))[i].L    = global::Ls + count*global::aLL;
	  (*Cells(T))[i].J    = global::Js + count*global::aPP;
	  count              += cpt;
        }else{
	  (*Cells(T))[i].L    = global::Ls + count*global::aLL;
	  (*Cells(T))[i].J    = global::Js + count*global::aPP;
	  count++;
	}
      }
    }
    for(int i = 0; i < internal::LsSize; i++){
      global::Ls[i] = cplx0;
    }
    for(int i = 0; i < internal::JsSize; i++){
      global::Js[i] = cplx0;
    }
  }



  
  template<int DIM>
  void InitializeDevs(Tree<DIM>& S, Tree<DIM>& T){
    // Allcoations for sources
    for(int i = 0; i < nCells(S); i++){
      if(global::precomputationAtStage[(*Cells(S))[i].lvl]){
	if(global::isHFStage[(*Cells(S))[i].lvl]){
	  int cpt = 0;
	  for(int k = 0; k < global::nDirPerStage[(*Cells(S))[i].lvl]; k++){
	    if((*Cells(S))[i].effDir[k]){
	      (*Cells(S))[i].effDir[k] = ++cpt;
	    }
	  }
	  (*Cells(S))[i].ndir = cpt;
	  (*Cells(S))[i].dirs = new int[(*Cells(S))[i].ndir];
	  int loccount = 0;
	  for(int u = 0; u < global::nDirPerStage[(*Cells(S))[i].lvl]; u++){
	    if((*Cells(S))[i].effDir[u]){(*Cells(S))[i].dirs[loccount++] = u;}
	  }
	  if(posix_memalign((void**)&((*Cells(S))[i].M), 64, sizeof(cplx)*cpt*global::aLL)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < cpt*global::aLL; j++){(*Cells(S))[i].M[j] = cplx0;}
	  if(posix_memalign((void**)&((*Cells(S))[i].F), 64, sizeof(cplx)*cpt*global::aPP)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < cpt*global::aPP; j++){(*Cells(S))[i].F[j] = cplx0;}
        }else{
	  if(posix_memalign((void**)&((*Cells(S))[i].M), 64, sizeof(cplx)*global::LL)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < global::LL; j++){(*Cells(S))[i].M[j] = cplx0;}
	  if(posix_memalign((void**)&((*Cells(S))[i].F), 64, sizeof(cplx)*global::PP)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < global::PP; j++){(*Cells(S))[i].F[j] = cplx0;}
	}
      }
    }
    // Allocate targets
    for(int i = 0; i < nCells(T); i++){
      if(global::precomputationAtStage[(*Cells(T))[i].lvl]){
	if(global::isHFStage[(*Cells(T))[i].lvl]){
	  int cpt = 0;
	  for(int k = 0; k < global::nDirPerStage[(*Cells(T))[i].lvl]; k++){
	    if((*Cells(T))[i].effDir[k]){
	      (*Cells(T))[i].effDir[k] = ++cpt;
	    }
	  }
	  (*Cells(T))[i].ndir = cpt;
	  (*Cells(T))[i].dirs = new int[(*Cells(T))[i].ndir];
	  int loccount = 0;
	  for(int u = 0; u < global::nDirPerStage[(*Cells(T))[i].lvl]; u++){
	    if((*Cells(T))[i].effDir[u]){(*Cells(T))[i].dirs[loccount++] = u;}
	  }
	  if(posix_memalign((void**)&((*Cells(T))[i].L), 64, sizeof(cplx)*cpt*global::aLL)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < cpt*global::aLL; j++){(*Cells(T))[i].L[j] = cplx0;}
	  if(posix_memalign((void**)&((*Cells(T))[i].J), 64, sizeof(cplx)*cpt*global::aPP)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < cpt*global::aPP; j++){(*Cells(T))[i].J[j] = cplx0;}
        }else{
	  if(posix_memalign((void**)&((*Cells(T))[i].L), 64, sizeof(cplx)*global::LL)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < global::LL; j++){(*Cells(T))[i].L[j] = cplx0;}
	  if(posix_memalign((void**)&((*Cells(T))[i].J), 64, sizeof(cplx)*global::PP)){std::cout << f_red << "Aligned alloc error: initdevs." << std::endl; exit(1);}	  
	  for(int j = 0; j < global::PP; j++){(*Cells(T))[i].J[j] = cplx0;}
	}
      }
    }
  }

  
} // DEFMM

#endif
