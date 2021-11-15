//===================================================================
//
// Authors: Igor Chollet, Xavier Claeys, Pierre Fortin, Laura Grigori
//
//  Copyright (C) 2020 Sorbonne Universite, Inria
//  All right reserved.
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

#ifndef INTERFACE_MAIN_HEADER_HPP
#define INTERFACE_MAIN_HEADER_HPP

#include <fstream>
#include "./structs/tree.hpp"
#include "./traversals/hrzpass.hpp"
#include "./traversals/upwpass.hpp"
#include "./traversals/dwnpass.hpp"
#include "./traversals/m2fpass.hpp"
#include "./traversals/f2lpass.hpp"
#include "./tables/precompute.hpp"
#include "./tables/frequency.hpp"
#include "./symmetries/permutations.hpp"
#include "global.hpp"
#include "./calculus/planewave.hpp"
#include "./tables/mapOnTrees.hpp"
#include "./structs/initdevs.hpp"

namespace defmm{

  template<int DIM>
  class IBFMM_Mat{
  private :
    int L;
    int Ncrit;
    flt kappa;
    Tree<DIM> S;
    Tree<DIM> T;
    particles<DIM> Y;
    particles<DIM> X;
    
  public :
    IBFMM_Mat(){}

    void addSourceParticle(flt part[DIM]){
      vecf<DIM> tmp;
      for(int d = 0; d < DIM; d++){
	tmp[d] = part[d];
      }
      particle<DIM> ttmp;
      ttmp.pos = tmp;
      Y.push_back(ttmp);
    }

    void addTargetParticle(flt part[DIM]){
      vecf<DIM> tmp;
      for(int d = 0; d < DIM; d++){
	tmp[d] = part[d];
      }
      particle<DIM> ttmp;
      ttmp.pos = tmp;
      X.push_back(ttmp);
    }

    void addSourceParticlesINP(const char* name, int nparts){
      std::ifstream infile(name);
      Y.resize(nparts);
      for(int i = 0; i < nparts; i++){
        for(int k = 0; k < DIM; k++){
	  infile >> Y[i].pos[k];
	}
      }
      infile.close();
    }
    void addTargetParticlesINP(const char* name, int nparts){
      std::ifstream infile(name);
      X.resize(nparts);
      for(int i = 0; i < nparts; i++){
	for(int k = 0; k < DIM; k++){
	  infile >> X[i].pos[k];
	}
      }
      infile.close();
    }

    void prcmpt(int _L, int _Ncrit, flt _kappa, int lvMAX = -1){
#ifdef DEFMM_VERBOSE
      flt t0,t1;
      t0 = gaia::chronos::time();
#endif
      L     = _L;
      Ncrit = _Ncrit;
      kappa = _kappa;
      T.setParticles(&X);
      S.setParticles(&Y);
      if(lvMAX == -1){
	DoubleBuild(T,S,L,Ncrit,kappa);
      }else{
	DoubleBuild(T,S,L,lvMAX);
      }
      global::q = new cplx[nParts(S)];
      global::p = new cplx[nParts(T)];
      global::setParameters<DIM>(L,kappa);
      computeZ2ZPermutationsForM2ML2L<DIM>();
      circulant::precomputeCirculant<DIM>();
      initFrequency(T,S);
      planewave::computeSignatures(T,S);
      allocatePermutations<DIM>();
      computeSigma<DIM>();
      computePermutations<DIM>();
      setHashchode2ActionsIndices<DIM>();
      allocateM2LTable<DIM,Vecc                      >(T,S);
      prcmputeM2LTable<DIM,thetis::kernels::Helmholtz>(T,S);
      fastlagrange::prcmptM2ML2L<DIM>();
      fastlagrange::realPrecomputeM2M<DIM>();
      fastlagrange::realPrecomputeL2L<DIM>();
      InitializeSourceArrays(S);
      InitializeTargetArrays(T);
      hrzpassinternal::cpt = 0;
      BlankHorizontalPass<DIM>(T,S);
      if(!global::noHF){
	BlankDownwardPass(S);
	BlankDownwardPass(T);
      }
#ifdef INITIALIZE_DEV_BLOCK
      InitializeSourceDevsBlock(S);
      InitializeTargetDevsBlock(T);
#else
      InitializeDevs(S,T);
#endif
      fastComputeSourceLeafMatrices(S);
      fastComputeTargetLeafMatrices(T);
#ifdef DEFMM_VERBOSE
      t1 = gaia::chronos::time();
      std::cout << "Precomputation time: " << f_cyan << t1-t0 << f_def << std::endl;
#endif
    }
    
    friend void gemv(IBFMM_Mat<DIM>& A, Vecc& In, Vecc& Out){
#ifdef DEFMM_VERBOSE
      flt t0,t1;
      t0 = gaia::chronos::time();
#endif
      A.S.loadRHS(In,global::q);
      for(int i = 0; i < nParts(A.T); i++){global::p[i] = cplx0;}
      UpwardPass(A.S);
      M2FPass<DIM>(A.S);
      HorizontalPass(A.T,A.S);
      F2LPass<DIM>(A.T);
      DownwardPass(A.T);
      A.T.writeLHS(Out,global::p);
#ifdef DEFMM_VERBOSE
      t1 = gaia::chronos::time();
      std::cout << "Application time: " << f_cyan << t1-t0 << f_def << std::endl;
#endif
    }
    
    void loadRHS(Vecc& In){S.loadRHS(In);}
    void writeLHS(Vecc& Out){T.writeLHS(Out);}
    friend cell<DIM>* SourceRoot(IBFMM_Mat& A){return &Root(A.S)[0];}
    friend cell<DIM>* TargetRoot(IBFMM_Mat& A){return &Root(A.T)[0];}


    friend void CheckFULL(IBFMM_Mat& A){
      Vecc Exa(nParts(A.T)); Exa = cplx0;
      std::cout << "Direct computation time: ";
      double t0 = gaia::chronos::time();
      int N = nParts(A.T);
      int M = nParts(A.S);
      flt* pt;
      flt* ps;
      flt  pi4m1 = 1./pi4;
      flt x, y, z, R, K, r, co, si, pr, pi, qr, qi;
      cplx *qq = global::q;
#pragma omp simd
      for(int i = 0; i < N; i++){
	pt = K_(Parts(A.T)[i].pos);
	pr = 0.;
	pi = 0.;
	for(int j = 0; j < M; j++){
	  ps     = K_(Parts(A.S)[j].pos);
	  qr     = qq[j].real();   // Deinterleave data
	  qi     = qq[j].imag();
	  x      = pt[0] - ps[0];
	  y      = pt[1] - ps[1];
	  z      = pt[2] - ps[2];
	  R      = x*x+y*y+z*z;
	  K      = 1./sqrt(R);     // Fast approx inverse sqrt
	  if(R < 1.e-16){K = 0.;}
	  r      = K * R;
	  r     *= global::kappa;
	  K     *= pi4m1;
	  co     = cos(r);
	  si     = sin(r);
	  co    *= K;
	  si    *= K;
	  pr    += co*qr - si*qi;
	  pi    += co*qi + si*qr;
	}
	Exa[i] += cplx(pr,pi);
      }
      double t1 = gaia::chronos::time();
      std::cout << f_cyan << t1 - t0 << f_def << std::endl;
      std::cout << "Relative L^oo norm: " << f_cyan << nrmInfRel(Exa,global::p) << f_def << std::endl;
      std::cout << "Relative L^1  norm: " << f_cyan << nrm1Rel  (Exa,global::p) << f_def << std::endl;
      std::cout << "Relative L^2  norm: " << f_cyan << nrm2Rel  (Exa,global::p) << f_def << std::endl;
    }

    
  };


} // DEFMM

#endif
