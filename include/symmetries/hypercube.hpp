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

#ifndef HYPERCUBE_HPP
#define HYPERCUBE_HPP

#include "../gaia/gaia.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{
  /*********************   Group represented as actions on R^d   *********************/
  template<int DIM>
  void HypercubeZo2Z(vecf<DIM>& v, int i){
    veci<DIM> I = thetis::keys::i2I<DIM>(i,2);
    for(int k = 0; k < DIM; k++){
      if(I[k]){v[k] *= -1;}
    }
  }
  template<int DIM>
  void HypercubeZo2Z(veci<DIM>& v, int i){
    veci<DIM> I = thetis::keys::i2I<DIM>(i,2);
    for(int k = 0; k < DIM; k++){
      if(I[k]){v[k] *= -1;}
    }
  }
  template<int DIM>
  void HypercubeZo2ZInvInd(int i, int& iz){
    iz = i;
  }
  template<int DIM>
  int HypercubeZo2ZCardinal(){
    return (1<<DIM);
  }
  template<int DIM>
  void HypercubeSigma(vecf<DIM>& v, int i){
    std::cout << f_red << "'HypercubeSigma' function (hypercube.hpp) is not implemented in the general case. You need to provide a specification for d = " << DIM << f_def << std::endl;
    exit(1);
  }
  template<int DIM>
  void HypercubeSigma(veci<DIM>& v, int i){
    std::cout << f_red << "'HypercubeSigma' function (hypercube.hpp) is not implemented in the general case. You need to provide a specification for d = " << DIM << f_def << std::endl;
    exit(1);
  }
  template<int DIM>
  void HypercubeSigmaInvInd(int i, int& is){
    std::cout << f_red << "'HypercubeSigmaInvInd' function (hypercube.hpp) is not implemented in the general case. You need to provide a specification for d = " << DIM << f_def << std::endl;
    exit(1);
  }
  template<int DIM>
  int HypercubeSigmaCardinal(){
    int p = 1;
    for(int i = 2; i <= DIM; i++){
      p *= i;
    }
    return p;
  }
  template<int DIM>
  int HypercubeCardinal(){
    return HypercubeZo2ZCardinal<DIM>()*HypercubeSigmaCardinal<DIM>();
  }

  template<int DIM>
  void HypercubeSymmetry(vecf<DIM>& v, int s, int z){
    HypercubeSigma(v,s);
    HypercubeZo2Z (v,z);
  }
  template<int DIM>
  void HypercubeSymmetry(veci<DIM>& v, int s, int z){
    HypercubeSigma(v,s);
    HypercubeZo2Z (v,z);
  }
  
  template<int DIM>
  void HypercubeSymmetry(vecf<DIM>& v, int i){
    std::cout << f_red << "'HypercubeSymmetry' function (permutation.hpp) is not implemented in the general case. You need to provide a specification for d = " << DIM << f_def << std::endl;
    //exit(1);
  }
  template<int DIM>
  void HypercubeSymmetry(veci<DIM>& v, int i){
    std::cout << f_red << "'HypercubeSymmetry' function (permutation.hpp) is not implemented in the general case. You need to provide a specification for d = " << DIM << f_def << std::endl;
    //exit(1);
  }
  
  /********************************************************************/
  // Hashcodes
  template<int DIM>
  inline int hashcodeSigma(veci<DIM>& I){
    int hashcode = 0, u = 0;
    for(int k = 0; k < DIM; k++){
      for(int l = k+1; l < DIM; l++){
        hashcode += ((std::abs(I[k]) < std::abs(I[l])) << u++);
      }
    }
    return hashcode;
  }
  template<int DIM>
  inline int hashcodeSigmaSingular(veci<DIM>& I){
    int hashcode = 0, u = 0;
    for(int k = 0; k < DIM; k++){
      for(int l = k+1; l < DIM; l++){
        hashcode += ((std::abs(I[k]) == std::abs(I[l])) << u++);
      }
    }
    return hashcode;
  }
  
  template<int DIM>
  inline int hashcodeZov2Z(veci<DIM>& I){
    int hashcode = 0;
    for(int k = 0; k < DIM; k++){
      hashcode += ((I[k] < 0) << k);
    }
    return hashcode;
  }
  
  template<int DIM>
  inline int hashcodeZov2ZSingular(veci<DIM>& I){
    int hashcode = 0;
    for(int k = 0; k < DIM; k++){
      hashcode += ((I[k] == 0) << k);
    }
    return hashcode;
  }
  template<int DIM>
  inline int hashcodeHypercube(veci<DIM>& I){
    return (hashcodeZov2Z(I) + (hashcodeSigma(I)<<DIM));
  }
  template<int DIM>
  inline int hashcodeSigma(vecf<DIM>& I){
    int hashcode = 0, u = 0;
    for(int k = 0; k < DIM; k++){
      for(int l = k+1; l < DIM; l++){
        hashcode += ((std::abs(I[k]) < std::abs(I[l])) << u++);
      }
    }
    return hashcode;
  }
  template<int DIM>
  inline int hashcodeZov2Z(vecf<DIM>& I){
    int hashcode = 0;
    for(int k = 0; k < DIM; k++){
      hashcode += ((I[k] < 0) << k);
    }
    return hashcode;
  }
  template<int DIM>
  inline int hashcodeHypercube(vecf<DIM>& I){
    return (hashcodeZov2Z(I) + (hashcodeSigma(I)<<DIM));
  }

  /********************************************************************/
  // Inversions
  template<int DIM>
  void HypercubeSymmetryInvInd(int s, int z, int& is, int& iz){
    HypercubeSigmaInvInd<DIM>(s,is);
    HypercubeZo2ZInvInd <DIM>(z,iz);
    veci<DIM> I = 1;
    HypercubeZo2Z (I,iz);
    HypercubeSigma(I,is);
    iz = hashcodeZov2Z(I);
  }

  
  /********************************************************************/
  // 1D case
  template<> void HypercubeSymmetry<1>(vecf<1>& v, int i){
    if(i == 0){return;}
    v[0] = -v[0];
  }
  template<> void HypercubeSymmetry<1>(veci<1>& v, int i){
    if(i == 0){return;}
    v[0] = -v[0];
  }
  template<> void HypercubeSigmaInvInd<1>(int i, int& is){
    is = i;
  }
  template<> void HypercubeSigma<1>(vecf<1>& v, int i){return;}
  template<> void HypercubeSigma<1>(veci<1>& v, int i){return;}
  template<> int HypercubeCardinal<1>(){return 2;}
  template<> int HypercubeZo2ZCardinal<1>(){return 2;}
  template<> int HypercubeSigmaCardinal<1>(){return 1;}

  /********************************************************************/
  // 2D case
  template<>
  void HypercubeSigma<2>(vecf<2>& v, int i){
    if(i == 0){return;}
    if(i == 1){flt tmp = v[0]; v[0] = v[1]; v[1] = tmp;}
  }
  template<>
  void HypercubeSigma<2>(veci<2>& v, int i){
    if(i == 0){return;}
    if(i == 1){flt tmp = v[0]; v[0] = v[1]; v[1] = tmp;}
  }
  template<> void HypercubeSigmaInvInd<2>(int i, int& is){
    is = i;
  }
  template<>
  int HypercubeCardinal<2>(){return 8;}
  template<> int HypercubeZo2ZCardinal<2>(){return 4;}
  template<> int HypercubeSigmaCardinal<2>(){return 2;}
  
  /********************************************************************/
  // 3D case
  template<>
  void HypercubeSigma<3>(vecf<3>& v, int i){
    if(i == 0){return;}
    if(i == 1){
      flt tmp = v[0];
      v[0] = v[1];
      v[1] = tmp;
      return;
    }
    if(i == 2){
      flt tmp = v[0];
      v[0] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 3){
      flt tmp = v[1];
      v[1] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 4){
      flt tmp = v[0];
      v[0] = v[1];
      v[1] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 5){
      flt tmp = v[0];
      v[0] = v[2];
      v[2] = v[1];
      v[1] = tmp;
      return;
    }
  }
  template<>
  void HypercubeSigma<3>(veci<3>& v, int i){
    if(i == 0){return;}
    if(i == 1){
      flt tmp = v[0];
      v[0] = v[1];
      v[1] = tmp;
      return;
    }
    if(i == 2){
      flt tmp = v[0];
      v[0] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 3){
      flt tmp = v[1];
      v[1] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 4){
      flt tmp = v[0];
      v[0] = v[1];
      v[1] = v[2];
      v[2] = tmp;
      return;
    }
    if(i == 5){
      flt tmp = v[0];
      v[0] = v[2];
      v[2] = v[1];
      v[1] = tmp;
      return;
    }
  }
  template<> void HypercubeSigmaInvInd<3>(int i, int& is){
    if(i == 0){is = 0; return;}
    if(i == 1){is = 1; return;}
    if(i == 2){is = 2; return;}
    if(i == 3){is = 3; return;}
    if(i == 4){is = 5; return;}
    if(i == 5){is = 4; return;}
  }
  template<>
  int HypercubeCardinal<3>(){return 48;}
  template<> int HypercubeZo2ZCardinal<3>(){return 8;}
  template<> int HypercubeSigmaCardinal<3>(){return 6;}
} // DEFMM

#endif
