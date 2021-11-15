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

#ifndef SORT_PARTICLES_HPP
#define SORT_PARTICLES_HPP

#include "particles.hpp"

#define STOP_NCRIT 1
#define STOP_LEVEL 2

namespace defmm{

  template<int DIM>
  inline void getBitcode(vecf<DIM>& u, vecf<DIM>& ref, int& bitcode){
    bitcode = 0;
    for(int k = 0; k < DIM; k++){
      bitcode += ((u[k] > ref[k] ? 1 : 0) << k);
    }
  }

  template<int STOP_ID>
  inline bool stopCriterion(int N, int lvl, int stopNumber){
    std::cout << f_red << "Warning: stopCriterion not defined for STOP_ID = " << STOP_ID << " in sort.hpp" << std::endl; exit(1);
  }
  template<>
  inline bool stopCriterion<STOP_NCRIT>(int N, int lvl, int stopNumber){return (N   < stopNumber);}
  template<>
  inline bool stopCriterion<STOP_LEVEL>(int N, int lvl, int stopNumber){return (lvl >= stopNumber);}
  
  /************************************************************************************/
  /******************************   Particle sort    **********************************/
  /************************************************************************************/
  
  /*
    parts : Vecteur des particules
    buff  : Buffer pour classer localement les particules
    pos   : Associe à une particule l'indice du 'fils' dans lequel elle se trouve
    N     : Nombre de particules dans le tableau local
    ctr   : Centre de la cellule locale
  */
  template<int DIM, int STOPCOND>
  void rcSort(particle<DIM>* parts, particle<DIM>* buff, int* pos, int N, vecf<DIM>& ctr, veci<DIM>& idx, flt rad, int lvl, int stopNumber, int& countCells){
    countCells++;
    if(stopCriterion<STOPCOND>(N,lvl,stopNumber)){return;}
    int  p2        = (1 << DIM);
    int* histo     = new int[p2]; // Va compter le nombre de particule dans chaque fils possible
    int* offsets   = new int[p2]; // Va stocker les offsets des fils
    int* addoffset = new int[p2];
    int  bitcode;
    for(int i = 0; i < p2; i++){histo[i] = 0; addoffset[i] = 0;}
    for(int i = 0; i < N; i++){
      getBitcode(parts[i].pos,ctr,bitcode);
      histo[bitcode]++;
      pos[i] = bitcode;
    }
    offsets[0] = 0;
    for(int i = 1; i < p2; i++){offsets[i] = histo[i-1]+offsets[i-1];}
    for(int i = 0; i < N; i++){buff[offsets[pos[i]]+(addoffset[pos[i]]++)] = parts[i];}
    for(int i = 0; i < N; i++){parts[i] = buff[i];}
    flt nextrad = rad/2.;
    for(int i = 0; i < p2; i++){
      if(histo[i]){
	vecf<DIM> nextctr = ctr;
	veci<DIM> nextidx = 0;
	for(int k = 0; k < DIM; k++){
	  nextctr[k] += nextrad * (((i & 1 << k) >> k) * 2 - 1);
	  nextidx[k]  = 2*idx[k] + ((i & (1<<k)) ? 1 : 0); // idx inutile...
        }
	rcSort<DIM,STOPCOND>(parts+offsets[i],buff+offsets[i],pos+offsets[i],histo[i],
			     nextctr,nextidx,nextrad,lvl+1,stopNumber,countCells);
      }
    }
  }

  template<int DIM, int STOPCOND>
  void rcSort(particle<DIM>* parts,
	      particle<DIM>* buff,
	      int* pos,
	      int N,
	      vecf<DIM>& ctr,
	      flt rad,
	      int lvl,
	      int stopNumber,
	      int& countCells){
    countCells++;
    // Stop the recursion?
    if(stopCriterion<STOPCOND>(N,lvl,stopNumber)){
      // Put the final result in the good array
      if((lvl%2)){
	for(int i = 0; i < N; i++){buff[i] = parts[i];}
      }
      return;
    }
    int  p2        = (1 << DIM);
    int* histo     = new int[p2]; // Number of particle in each son
    int* offsets   = new int[p2]; // Store son's offsets
    int* addoffset = new int[p2]; // Cumulative offsets
    int  bitcode;
    for(int i = 0; i < p2; i++){histo[i] = 0; addoffset[i] = 0;}
    // Decide for each particle in which son it is
    for(int i = 0; i < N; i++){
      bitcode = 0;
      for(int k = 0; k < DIM; k++){
	bitcode += ((parts[i].pos[k] > ctr[k]) << k);
      }
      histo[bitcode]++;
      pos[i] = bitcode;
    }
    // Set offsets
    offsets[0] = 0;
    for(int i = 1; i < p2; i++){
      offsets[i] = histo[i-1]+offsets[i-1];
    }
    // Move the particles to the particle buffer
    // according to their offset
    for(int i = 0; i < N ; i++){
      buff[offsets[pos[i]]+(addoffset[pos[i]]++)] = parts[i];
    }
    // Proceed recursively for each non-empty son
    flt nextrad = rad/2.;
    for(int i = 0; i < p2; i++){
      if(histo[i]){
	vecf<DIM> nextctr = ctr;
	for(int k = 0; k < DIM; k++){
	  nextctr[k] += nextrad * (((i & 1 << k) >> k) * 2 - 1);
        }
	// Interleave buff and part argument position !!
	rcSort<DIM,STOPCOND>(buff+offsets[i],
			     parts+offsets[i],
			     pos+offsets[i],
			     histo[i],
			     nextctr,
			     nextrad,
			     lvl+1,
			     stopNumber,
			     countCells);
      }
    }
  }
  
  template<int DIM, int STOPCOND>
  void sort(particle<DIM>* parts, int N, vecf<DIM>& boxCtr, flt boxRad, int stopNumber, int& countCells){
    particles<DIM>   buff(N);
    std::vector<int> pos(N);
    veci<DIM>        idx = 0;
    rcSort<DIM, STOPCOND>(parts, &buff[0], &pos[0], N, boxCtr, boxRad, 0, stopNumber, countCells);
  }

  /************************************************************************************/
  /****************************   Cell construction  **********************************/
  /************************************************************************************/
  
  template<int DIM, template<int CELLTYPE_DIM> class CELLTYPE>
  inline void createCell(CELLTYPE<DIM>& c, veci<DIM>& idx, int lvl, int nprt, vecf<DIM>& ctr, flt& rad, particle<DIM>* prt, int ind){
    c.idx  = idx ;
    c.lvl  = lvl ;
    c.nprt = nprt;
    c.ctr  = ctr ;
    c.rad  = rad ;
    c.prt  = prt ;
    c.ind  = ind;
  }
  
  template<int DIM, int STOPCOND, template<int CELLTYPE_DIM> class CELLTYPE>
  void rcCreateCells(CELLTYPE<DIM>* cellarray, particle<DIM>* parts, int N, vecf<DIM>& boxCtr, flt& boxRad,veci<DIM>& idx, int lvl, int stopNumber, int& globalInd, int particleInd, std::vector<int>& nCellPerStage){
    int locInd = globalInd;
    createCell<DIM,CELLTYPE>(cellarray[locInd],idx,lvl,N,boxCtr,boxRad,parts,particleInd);
    if(nCellPerStage.size() <= (unsigned int)lvl){
      nCellPerStage.push_back(0);}
    nCellPerStage[lvl]++;
    if(stopCriterion<STOPCOND>(N,lvl,stopNumber)){return;}
    int  p2         = (1<<DIM);
    int* histo      = new int[p2];
    int* offsets    = new int[p2];
    int  bitcode;
    for(int i = 0; i < p2; i++){histo[i] = 0;}
    for(int i = 0; i < N; i++){
      getBitcode(parts[i].pos,boxCtr,bitcode);
      histo[bitcode]++;
    }
    offsets[0] = 0;
    for(int i = 1; i < p2; i++){offsets[i] = histo[i-1]+offsets[i-1];}
    flt nextrad = boxRad/2.;
    for(int i = 0; i < p2; i++){
      if(histo[i]){
        vecf<DIM> nextctr = boxCtr;
	veci<DIM> nextidx;
	for(int k = 0; k < DIM; k++){
	  nextctr[k] += nextrad * (((i & 1 << k) >> k) * 2 - 1);
	  nextidx[k]  = 2*idx[k] + ((i & (1<<k)) ? 1 : 0);
	}
	cellarray[locInd].son.push_back(++globalInd);
	rcCreateCells<DIM,STOPCOND,CELLTYPE>(cellarray,parts+offsets[i],histo[i],
					     nextctr,nextrad,nextidx,lvl+1,stopNumber,globalInd,particleInd+offsets[i], nCellPerStage);
      }
    }
  }

  template<int DIM, int STOPCOND, template<int CELLTYPE_DIM> class CELLTYPE>
  void createCells(std::vector<CELLTYPE<DIM> >& cellarray, int nCells, particle<DIM>* parts, int N, vecf<DIM>& boxCtr, flt boxRad, int stopNumber, std::vector<int>& nCellPerStage){
    cellarray.resize(nCells);
    veci<DIM> id0  = 0;
    int       ind0 = 0;
    rcCreateCells<DIM,STOPCOND,CELLTYPE>(&cellarray[0], parts, N, boxCtr, boxRad, id0, 0, stopNumber,ind0,0, nCellPerStage);
  }  
} // DEFMM


#endif
