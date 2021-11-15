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

#ifndef DIRECTIONS_HEADER_HPP
#define DIRECTIONS_HEADER_HPP

#include <bitset>
#include "../symmetries/hypercube.hpp"


namespace defmm{

  const flt OneMinusEpsilon = 1.-1.e-12;
  
  /*
    Get directions from refinement of the fundamental domain of \frak{D} in any dimension
   */
  template<int DIM>
  struct BoundarySimplex{
    vecf<DIM> vertices[DIM];
  }; // BoundarySimplex

  // Number of sons for a simplex refinement
  const int nSonsBoundarySimplex[4] = {-1,1,2,4};
  
  template<int DIM>
  struct DirectionTree{
    vecf<DIM>             direction;
    int                   index;
  };

  template<int DIM>
  struct DirectionTables{
    DirectionTree<DIM>** tables;
    int*                 indices;
    int*                 nsons;
  };
  
  ///////////////////////////////////////////////////////
  
  template<int DIM>
  void initDirectionTables(DirectionTables<DIM>& DirTab, int nLvl){
    DirTab.tables  = new DirectionTree<DIM>*[nLvl];
    DirTab.indices = new int[nLvl];
    DirTab.nsons   = new int[nLvl];
    for(int i = 0; i < nLvl; i++){
      DirTab.indices[i] = 0;
      int npi = (1<<((DIM-1)*i))*2*DIM;
      DirTab.tables[i] = new DirectionTree<DIM>[npi];
    }
  }
  
  template<int DIM>
  void iterativeGetDirections(DirectionTables<DIM>& Directions, int nStages){

    // Allocate tables
    initDirectionTables(Directions,nStages);

    // Get first reference directions
    for(int k = 0; k < DIM; k++){
      Directions.tables[0][2*k  ].index        =  2*k;
      Directions.tables[0][2*k  ].direction    =  0.;
      Directions.tables[0][2*k  ].direction[k] =  1.; // d = ( 1, 0, ..., 0)
      Directions.tables[0][2*k+1].index        =  2*k+1;
      Directions.tables[0][2*k+1].direction    =  0.;
      Directions.tables[0][2*k+1].direction[k] = -1.; // d = (-1, 0, ..., 0)
    }

    // Get hyperplan definition for faces
    vecf<DIM>* hypfaces = new vecf<DIM>[(1<<(DIM-1))];
    for(int k = 0; k < (1<<(DIM-1)); k++){
      hypfaces[k]    =  1.;
      hypfaces[k][0] =  0.;
      for(int q = 0; q < DIM-1; q++){
	if((k&(1<<q))){
	  hypfaces[k][1+q] = -1.;
	}
      }
    }

    // Number of sons per direction
    int nson = (1<<(DIM-1));

    // Get directions in a face of the hypercube
    int nson_pow_em1 = 1;
    for(int e = 1; e < nStages; e++){
      flt radius = 1./flt((1<<(e+1)));
      for(int i = 0; i < nson_pow_em1; i++){
	for(int q = 0; q < nson; q++){
	  Directions.tables[e][i*nson+q].direction = Directions.tables[e-1][i].direction;
	  for(int k = 0; k < DIM; k++){
	    Directions.tables[e][i*nson+q].direction[k] += hypfaces[q][k] * radius;
	  }
	  Directions.tables[e][i*nson+q].direction /= norm2(Directions.tables[e][i*nson+q].direction);
        }
      }
      nson_pow_em1 *= nson;
    }

    // Get direction on the opposit face using reflection
    int nson_pow_e = 1;
    for(int e = 1; e < nStages; e++){
      nson_pow_e *= nson;
      int faceSize = nson_pow_e;
      for(int i = 0; i < faceSize; i++){
	for(int q = 0; q < nson; q++){
	  Directions.tables[e][faceSize+i*nson+q].direction = Directions.tables[e][i*nson+q].direction;
	  Directions.tables[e][faceSize+i*nson+q].direction[0] *= -1.;
	}
      }
    }

    // Get other directions using rotations
    nson_pow_e = 1;
    for(int e = 1; e < nStages; e++){
      nson_pow_e *= nson;
      int faceSize2 = 2*nson_pow_e;
      for(int j = 0; j < DIM; j++){
	for(int ii = 0; ii < faceSize2; ii++){
	  for(int k = 0; k < DIM; k++){
	    int l = (((k+j)%DIM)>=0 ? ((k+j)%DIM) : ((k+j)%DIM)+DIM);
	    Directions.tables[e][j*faceSize2+ii].direction[k] = Directions.tables[e][ii].direction[l];
	  }
	}
      }
    }

    // Set indices component
    nson_pow_e = 1;
    for(int e = 1; e < nStages; e++){
      nson_pow_e *= nson;
      int faceSizeTot = 2*DIM*nson_pow_e;
      int count = 0;
      for(int ii = 0; ii < faceSizeTot; ii++){
	Directions.tables[e][ii].index = count++;
      }
    }
    
    // SET "INDICES" ATTRIBUTE:
    nson_pow_e = 1;
    for(int e = 0; e < nStages; e++){
      Directions.indices[e] = nson_pow_e*2*DIM;
      nson_pow_e *= nson;
    }
    // SET "NSONS" ATTRIBUTE:
    for(int i = 0; i < nStages; i++){
      Directions.nsons[i] = nson;
    }
  }

  
  template<int DIM>
  int bestConeIndex(DirectionTables<DIM>& dt, int nStage, vecf<DIM>& dir){
    flt tmp = 1.e30;
    int ind = -1;
    for(int i = 0; i < dt.indices[nStage]; i++){
      if(tmp > norm2(dt.tables[nStage][i].direction - dir)){
	tmp = norm2(dt.tables[nStage][i].direction - dir);
	ind = i;
      }
    }
    return ind;
  }

  template<int DIM>
  int fastBestConeIndex(DirectionTables<DIM>& dt, int nStage, vecf<DIM>& dir){
    flt tmp  = 1.e30;
    int ind;
    int nson = (1 << (DIM-1));
    int curind;
    curind = 0;
    for(int k = 0; k <= nStage; k++){
      tmp = 1.e30;
      for(int i = 0; i < nson; i++){
	if(tmp > norm2(dt.tables[k][curind+i].direction - dir)){
	  tmp = norm2(dt.tables[k][curind+i].direction - dir);
	  ind = i;
	}
      }
      curind = (curind+ind)*nson;
    }
    return ind;
  }
  
} // DEFMM

#endif
