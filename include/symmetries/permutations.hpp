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

#ifndef PERMUTATIONS_MAIN_HEADER_HPP
#define PERMUTATIONS_MAIN_HEADER_HPP

#include <bits/stdc++.h>
#include <vector>
#include "../gaia/gaia.hpp"
#include "../global.hpp"
#include "hypercube.hpp"

namespace defmm{

  /*
    This file contains permutations 
    for action of the hypercube group
    in any dimension.
   */

  namespace permutations{
    std::vector<int*> sigma;
    std::vector<Veci> table;
    std::vector<int*> zo2z_L;
    int               identity;
    int               nPerm;
    int*              hashcode2action;
  }

  template<int DIM>
  void allocatePermutations(){
    permutations::nPerm = (1<<(DIM+DIM*(DIM-1)/2));
    permutations::hashcode2action = new int[permutations::nPerm];
    permutations::table.resize(permutations::nPerm);
    permutations::sigma.resize((1<<(DIM*(DIM-1)/2)));
    for(int i = 0; i < permutations::nPerm; i++){
      permutations::table[i].alignAllocate(global::PP,64);
    }
    for(int i = 0; i < (1<<(DIM*(DIM-1)/2)); i++){
      permutations::sigma[i] = new int[DIM];
    }
  }

  template<int DIM>
  void computeSigma(){
    std::cout << "'definePermutations' needs to be specified for DIM = "
	      << DIM
	      << ". Exit(1)."
	      << std::endl;
    exit(1);
  }

  template<>
  void computeSigma<1>(){permutations::sigma[0][0] = 0;}

  template<>
  void computeSigma<2>(){
    permutations::sigma[0][0] = 0; // i >= j
    permutations::sigma[0][1] = 1;
    //--------------------
    permutations::sigma[1][0] = 1; // i < j
    permutations::sigma[1][1] = 0;
  }
  
  template<>
  void computeSigma<3>(){
    permutations::sigma[0][0] = 0; // i=j=k -> Identity (always)
    permutations::sigma[0][1] = 1;
    permutations::sigma[0][2] = 2;
    //--------------------
    permutations::sigma[1][0] = 1; // i<j; i>=k; j>=k
    permutations::sigma[1][1] = 0;
    permutations::sigma[1][2] = 2;
    //--------------------
    permutations::sigma[2][0] = 0; // Does not exist
    permutations::sigma[2][1] = 1;
    permutations::sigma[2][2] = 2;
    //--------------------
    permutations::sigma[3][0] = 1; // i<j; i<k; j>=k
    permutations::sigma[3][1] = 2;
    permutations::sigma[3][2] = 0;
    //--------------------
    permutations::sigma[4][0] = 0; // i>=j; i>=k; j<k
    permutations::sigma[4][1] = 2;
    permutations::sigma[4][2] = 1;
    //--------------------
    permutations::sigma[5][0] = 0; // Does not exist
    permutations::sigma[5][1] = 1;
    permutations::sigma[5][2] = 2;
    //--------------------
    permutations::sigma[6][0] = 2; // i>=j; i<k; j<k
    permutations::sigma[6][1] = 0;
    permutations::sigma[6][2] = 1;
    //--------------------
    permutations::sigma[7][0] = 2; // i<j; i<k; j<k
    permutations::sigma[7][1] = 1;
    permutations::sigma[7][2] = 0;
  }

  /************** M2L permutations **************/
  // These permutations act on the full hypercube
  // group.
  
  // Apply sign permut' to sigma's one.
  template<int DIM>
  void computePermutations(){
    for(int p0 = 0; p0 < (1 << DIM); p0++){
      for(int p1 = 0; p1 < (1 << (DIM*(DIM-1)/2)); p1++){
	int hashcode = p0 + (p1 << DIM);
        // Compute efective permutation, combining
	// sigma and Z/2Z...
	for(int p = 0; p < global::PP; p++){
	  veci<DIM> tp = thetis::keys::i2I<DIM>(p,global::P);
	  veci<DIM> tm;
	  for(int k = 0; k < DIM; k++){
	    tm[k] = tp[permutations::sigma[p1][k]];
	    if(p0 & (1<<permutations::sigma[p1][k])){
	      if(tm[k]){
		tm[k] = global::P - tm[k];
	      }
	    }
	  }
	  permutations::table[hashcode][p] = thetis::keys::I2i<DIM>(tm,global::P);
	}
      }
    }
    permutations::identity = 0;
  }

  /************** M/L permutations **************/
  // Only acting on the abelian subgroup of the
  // hypercube one : (Z/2Z)^d
  // => only the corresponding hashcode is needed
  template<int DIM>
  void computeZ2ZPermutationsForM2ML2L(){
    int nZperms = (1<<DIM);
    permutations::zo2z_L.resize(nZperms);
    for(int i = 0; i < nZperms; i++){
      permutations::zo2z_L[i] = new int[global::LL];
      for(int l = 0; l < global::LL; l++){
	veci<DIM> tl = thetis::keys::i2I<DIM>(l,global::L);
        for(int k = 0; k < DIM; k++){
	  if(i & (1<<k)){
	    tl[k] = global::L - tl[k] - 1;
	  }
	}
	permutations::zo2z_L[i][l] = thetis::keys::I2i<DIM>(tl,global::L);
      }
    }
  }

  template<int DIM>
  void setHashchode2ActionsIndices(){
    vecf<DIM> t;
    for(int k = 0; k < DIM; k++){t[k] = flt(DIM-k);}
    for(int i = 0; i < permutations::nPerm; i++){
      permutations::hashcode2action[i] = -1;
    }
    for(int z = 0; z < HypercubeZo2ZCardinal<DIM>(); z++){
      for(int s = 0; s < HypercubeSigmaCardinal<DIM>(); s++){
	vecf<DIM> tmp = t;
	HypercubeSymmetry<DIM>(tmp,s,z);
	int hashcode = hashcodeZov2Z<DIM>(tmp) + (hashcodeSigma<DIM>(tmp) << DIM);
	permutations::hashcode2action[hashcode] = z*HypercubeSigmaCardinal<DIM>()+s;
      }
    }
  }

  // Returns the permutation corresponding to
  // a given translation vector.
  template<int DIM>
  inline int* getPermute(vecf<DIM>& T){
    int hashcode = hashcodeZov2Z(T) + (hashcodeSigma(T) << DIM);
    return K_(permutations::table[hashcode]);
  }

  template<int DIM>
  inline int getPermuteHsh(veci<DIM>& T){
    int hashcode = hashcodeZov2Z(T) + (hashcodeSigma(T) << DIM);
    return hashcode;
  }

  // Returns the M2L vector of the multi-id
  // before permutation, i.e. the M2L
  // vector in the fundamental domain.
  template<int DIM>
  inline void getPermuteAndFundamentalMultiIndex(veci<DIM>& X, veci<DIM>& R, int*& permut){
    int hashcode0 = hashcodeZov2Z(X);
    int hashcode1 = hashcodeSigma(X);
    int hashcode2 = hashcode0 + (hashcode1<<DIM);
    // First Sigma and then Z/2Z
    int* per = permutations::sigma[hashcode1];
    for(int k = 0; k < DIM; k++){R[k] = std::abs(X[per[k]]);}
    // Return pointer on the corresponding permutation
    permut = K_(permutations::table[hashcode2]);
  }
  
  // Aslo returns the hashcode
  template<int DIM>
  inline void getPermuteAndFundamentalMultiIndexAndHashcode(veci<DIM>& X, veci<DIM>& R, int*& permut, int& hsh){
    int hashcode0 = hashcodeZov2Z(X);
    int hashcode1 = hashcodeSigma(X);
    int hashcode2 = hashcode0 + (hashcode1<<DIM);
    // First Sigma and then Z/2Z
    int* per = permutations::sigma[hashcode1];
    for(int k = 0; k < DIM; k++){R[k] = std::abs(X[per[k]]);}
    // Return pointer on the corresponding permutation
    permut = K_(permutations::table[hashcode2]);
    hsh = hashcode2;
  }

} // DEFMM


#endif
