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

#ifndef TREES_HEADER_HPP
#define TREES_HEADER_HPP

#include <stdint.h>
#include <cassert>
#include "particles.hpp"
#include "sort.hpp"

namespace defmm{

  template<int DIM>
  class Tree{
  private :
    vecf<DIM>         BoxCenter;
    flt               BoxRadius;
    cells<DIM>        tree;
    particles<DIM>*   tmpparts;
    particles<DIM>    parts;
    int               Ncrit;
    std::vector<int>  nCellPerStage;
    flt               leafRad;
  public  :
    Tree(){}
    Tree(particles<DIM>* _tmpparts) : tmpparts(_tmpparts) {};
    ~Tree(){
      parts.clear();
      tree.clear();
    }
    void setParticles(particles<DIM>* _tmpparts){tmpparts = _tmpparts;}

    void setLeafRad(flt lr){leafRad = lr;}

    inline friend cell<DIM>*      Root  (Tree<DIM>& T){return &T.tree[0];}
    inline friend int             nCells(Tree<DIM>& T){return T.tree.size();}
    inline friend int             nParts(Tree<DIM>& T){return T.parts.size();}
    inline friend particle<DIM>*  Parts (Tree<DIM>& T){return &T.parts[0];}
    inline friend int             Nlvl  (Tree<DIM>& T){return T.nCellPerStage.size();}
    inline friend flt             Radius(Tree<DIM>& T){return T.BoxRadius;}
    inline friend cells<DIM>*     Cells (Tree<DIM>& T){return &T.tree;}
    
    friend void SetParticleIdx(Tree<DIM>& T){
      for(unsigned int i = 0; i < (*T.tmpparts).size(); i++){
	(*T.tmpparts)[i].ind = i;
      }
    }

    friend void DoubleBuild(Tree<DIM>& T, Tree<DIM>& S, int L, int _Ncrit, flt kappa){
      SetParticleIdx(T);
      SetParticleIdx(S);
      T.Ncrit = _Ncrit;
      S.Ncrit = _Ncrit;
      T.parts.resize(T.tmpparts->size());
      T.tree.resize(0);
      S.parts.resize(S.tmpparts->size());
      S.tree.resize(0);
      getBounds(*T.tmpparts,T.BoxCenter,T.BoxRadius,*S.tmpparts,S.BoxCenter,S.BoxRadius);
      veci<DIM> id0 = 0;
      T.leafRad = T.BoxRadius;
      S.leafRad = S.BoxRadius;
      int cptcells = 0;
      for(unsigned int i = 0; i < T.tmpparts->size(); i++){
	T.parts[i] = (*T.tmpparts)[i];
      }
      if(!(*T.tmpparts).size()){std::cout << "Warning: N = 0" << std::endl;}
      sort<DIM,STOP_NCRIT>(&(T.parts)[0],T.parts.size(),T.BoxCenter,T.BoxRadius,_Ncrit,cptcells);
      createCells<DIM,STOP_NCRIT,cell>(T.tree, cptcells, &T.parts[0], T.parts.size(),T.BoxCenter,T.BoxRadius,_Ncrit,T.nCellPerStage);
      cptcells = 0;
      for(unsigned int i = 0; i < S.tmpparts->size(); i++){
	S.parts[i] = (*S.tmpparts)[i];
      }
      if(!(*S.tmpparts).size()){std::cout << "Warning: N = 0" << std::endl;}
      sort<DIM,STOP_NCRIT>(&(S.parts)[0],S.parts.size(),S.BoxCenter,S.BoxRadius,_Ncrit,cptcells);
      createCells<DIM,STOP_NCRIT,cell>(S.tree, cptcells, &S.parts[0], S.parts.size(),S.BoxCenter,S.BoxRadius,_Ncrit,S.nCellPerStage);
      (*T.tmpparts).clear();
      id0 = 0;
      (*S.tmpparts).clear();
    }
    
    friend void DoubleBuild(Tree<DIM>& T, Tree<DIM>& S, int L, int _lvMax){
      SetParticleIdx(T);
      SetParticleIdx(S);
      T.parts.resize(T.tmpparts->size());
      T.tree.resize(0);
      S.parts.resize(S.tmpparts->size());
      S.tree.resize(0);
      getBounds(*T.tmpparts,T.BoxCenter,T.BoxRadius,*S.tmpparts,S.BoxCenter,S.BoxRadius);
      veci<DIM> id0 = 0;
      T.leafRad = T.BoxRadius;
      S.leafRad = S.BoxRadius;
      int cptcells = 0;
      for(unsigned int i = 0; i < T.tmpparts->size(); i++){
	T.parts[i] = (*T.tmpparts)[i];
      }
      if(!(*T.tmpparts).size()){std::cout << "Warning: N = 0" << std::endl;}
      sort<DIM,STOP_LEVEL>(&(T.parts)[0],T.parts.size(),T.BoxCenter,T.BoxRadius,_lvMax,cptcells);
      createCells<DIM,STOP_LEVEL,cell>(T.tree, cptcells, &T.parts[0], T.parts.size(),T.BoxCenter,T.BoxRadius,_lvMax,T.nCellPerStage);
      cptcells = 0;
      for(unsigned int i = 0; i < S.tmpparts->size(); i++){
	S.parts[i] = (*S.tmpparts)[i];
      }
      if(!(*S.tmpparts).size()){std::cout << "Warning: N = 0" << std::endl;}
      sort<DIM,STOP_LEVEL>(&(S.parts)[0],S.parts.size(),S.BoxCenter,S.BoxRadius,_lvMax,cptcells);
      createCells<DIM,STOP_LEVEL,cell>(S.tree, cptcells, &S.parts[0], S.parts.size(),S.BoxCenter,S.BoxRadius,_lvMax,S.nCellPerStage);
      (*T.tmpparts).clear();
      id0 = 0;
      (*S.tmpparts).clear();
    }

    friend void InitializeDevs(Tree<DIM>& T, int L){
      int l_pow_dim = 1;
      int p_pow_dim = 1;
      for(int i = 1; i <= DIM; i++){
	l_pow_dim *= L;
	p_pow_dim *= 2*L-1;
      }
      for(int i = 0; i < T.tree.size(); i++){
	T.tree[i].M.allocate(l_pow_dim);     T.tree[i].M = cplx0;
	T.tree[i].L.allocate(l_pow_dim);     T.tree[i].L = cplx0;
	T.tree[i].F.allocate(p_pow_dim); T.tree[i].F = cplx0;
	T.tree[i].J.allocate(p_pow_dim); T.tree[i].J = cplx0;
      }
    }

    inline void loadRHS(Vecc& Q, cplx* q){
      for(unsigned int i = 0; i < parts.size(); i++){
	q[i] = Q[parts[i].ind];
      }

    }

    inline void writeLHS(Vecc& P, cplx* p){
      for(unsigned int i = 0; i < parts.size(); i++){
	P[parts[i].ind] = p[i];
      }
    }

    friend void printTreeInfos(Tree<DIM>& T){
      std::cout << "Tree infos :\t(Nb Stages) " << f_cyan << T.nCellPerStage.size() << f_def
		<< ";\t(Nb Cells) " << f_cyan << T.tree.size() << f_def << std::endl;
    }
    
  };
  
} // DEFMM

#endif
