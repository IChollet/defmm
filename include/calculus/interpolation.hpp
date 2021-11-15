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

#ifndef DUFMM_INTERPOLATION_MAIN_HEADER_HPP
#define DUFMM_INTERPOLATION_MAIN_HEADER_HPP

#include "../structs/particles.hpp"
#include "../global.hpp"
#include "../thetis/thetis.hpp"

namespace defmm{

  namespace fastlagrange{

    /*
      Buffer stuff for single expansion
    */
    inline void loadBuf(cplx* M_s_dev, Vecc& bufq){
      for(int i = 0; i < global::LL; i++){
	bufq[i] = M_s_dev[i];
      }
    }
    inline void writeBuf(cplx* M_t_dev, Vecc& bufq){
      for(int i = 0; i < global::LL; i++){
	M_t_dev[i] += bufq[i];
      }
    }
    inline void loadBuf(cplx* M_s_dev, cplx* bufq){
      for(int i = 0; i < global::LL; i++){
	bufq[i] = M_s_dev[i];
      }
    }
    inline void writeBuf(cplx* M_t_dev, cplx* bufq){
      for(int i = 0; i < global::LL; i++){
	M_t_dev[i] += bufq[i];
      }
    }

    namespace internal{
      Matc* M2M_S1DxS1D;
      Matf* real_M2M_S1DxS1D;
      int*  M2M_permutation;
      Matc* sparse_M2M_S1DxS1D;
      int*  sparse_M2M_permutation_trg;
      int*  sparse_M2M_permutation_src;
      int   sparse_M2M_size_trg;
      int   sparse_M2M_size_src;
      Matc* L2L_S1DxS1D;
      Matf* real_L2L_S1DxS1D;
      int*  L2L_permutation;
      int   LpDm1;
    } // INTERNAL

    template<int DIM>
    void precomputeM2M(){
      internal::M2M_S1DxS1D = new Matc[2];
      internal::M2M_S1DxS1D[0].allocate(global::L,global::L);
      internal::M2M_S1DxS1D[1].allocate(global::L,global::L);
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::M2M_S1DxS1D[0](k,p) = thetis::uniform::S1D(k,
							       thetis::uniform::point1D(p,
											global::L,
											1.)-1.,
							       global::L,
							       2.);
	} 
      }
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::M2M_S1DxS1D[1](k,p) = thetis::uniform::S1D(k,
							       thetis::uniform::point1D(p,
											global::L,
											1.)+1.,
							       global::L,
							       2.);
	} 
      }
      internal::M2M_permutation = new int[global::LL];
      for(int i = 0; i < global::LL; i++){
	veci<DIM> I  = thetis::keys::i2I<DIM>(i,global::L);
	int     tmp  = I[0];
	for(int k = 1; k < DIM; k++){
	  I[k-1] = I[k];
	}
	I[DIM-1] = tmp;
	int  itilda  = thetis::keys::I2i<DIM>(I,global::L);
	internal::M2M_permutation[i] = itilda;
      }
    }

    template<int DIM>
    void precomputeL2L(){
      internal::L2L_S1DxS1D = new Matc[2];
      internal::L2L_S1DxS1D[0].allocate(global::L,global::L);
      internal::L2L_S1DxS1D[1].allocate(global::L,global::L);
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::L2L_S1DxS1D[0](p,k) = thetis::uniform::S1D(k,
							       thetis::uniform::point1D(p,
											global::L,
											1.)-1.,
							       global::L,
							       2.);
	} 
      }
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::L2L_S1DxS1D[1](p,k) = thetis::uniform::S1D(k,
							       thetis::uniform::point1D(p,
											global::L,
											1.)+1.,
							       global::L,
							       2.);
	} 
      }
      internal::L2L_permutation = new int[global::LL];
      for(int i = 0; i < global::LL; i++){
	veci<DIM> I  = thetis::keys::i2I<DIM>(i,global::L);
	int     tmp  = I[0];
	for(int k = 1; k < DIM; k++){
	  I[k-1] = I[k];
	}
	I[DIM-1] = tmp;
	int  itilda  = thetis::keys::I2i<DIM>(I,global::L);
	internal::L2L_permutation[i] = itilda;
      }
    }

    template<int DIM>
    void precomputeSparseM2M(){
      internal::sparse_M2M_size_src = global::L/2;
      internal::sparse_M2M_size_trg = global::L - internal::sparse_M2M_size_src;
      internal::sparse_M2M_permutation_src = new int[internal::sparse_M2M_size_src];
      internal::sparse_M2M_permutation_trg = new int[internal::sparse_M2M_size_trg];
      internal::sparse_M2M_S1DxS1D = new Matc[2];
      internal::sparse_M2M_S1DxS1D[0].allocate(global::L,internal::sparse_M2M_size_src);
      internal::sparse_M2M_S1DxS1D[1].allocate(global::L,internal::sparse_M2M_size_src);
      // Get offsets for a vector of size L
      int countsrc = 0, counttrg = 0;
      for(int i = 0; i < global::L; i++){
	if(i%2){
	  internal::sparse_M2M_permutation_src[countsrc++] = i;
	}else{
	  internal::sparse_M2M_permutation_trg[counttrg++] = i;
	}
      }
      // Get sparse M2M matrices
      for(int p = 1; p < global::L; p+=2){
	for(int k = 0; k < global::L; k++){
	  internal::sparse_M2M_S1DxS1D[0](k,p/2) = thetis::uniform::S1D(k,thetis::uniform::point1D(p,global::L,1.)-1.,global::L,2.);
	  if((global::L%2)){
	    internal::sparse_M2M_S1DxS1D[1](k,p/2) = thetis::uniform::S1D(k,thetis::uniform::point1D(p,global::L,1.)+1.,global::L,2.);
	  }else{
	    internal::sparse_M2M_S1DxS1D[1](k,p/2) = thetis::uniform::S1D(k,thetis::uniform::point1D(p-1,global::L,1.)+1.,global::L,2.);
	  }
	} 
      }
    }

    template<int DIM>
    void realPrecomputeM2M(){
      internal::real_M2M_S1DxS1D = new Matf[2];
      internal::real_M2M_S1DxS1D[0].allocate(global::L,global::L);
      internal::real_M2M_S1DxS1D[1].allocate(global::L,global::L);
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::real_M2M_S1DxS1D[0](k,p) = thetis::uniform::S1D(k,
								    thetis::uniform::point1D(p,
											     global::L,
											     1.)-1.,
								    global::L,
								    2.);
	} 
      }
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::real_M2M_S1DxS1D[1](k,p) = thetis::uniform::S1D(k,
								    thetis::uniform::point1D(p,
											     global::L,
											     1.)+1.,
								    global::L,
								    2.);
	} 
      }
    }

    template<int DIM>
    void realPrecomputeL2L(){
      internal::real_L2L_S1DxS1D = new Matf[2];
      internal::real_L2L_S1DxS1D[0].allocate(global::L,global::L);
      internal::real_L2L_S1DxS1D[1].allocate(global::L,global::L);
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::real_L2L_S1DxS1D[0](p,k) = thetis::uniform::S1D(k,
								    thetis::uniform::point1D(p,
											     global::L,
											     1.)-1.,
								    global::L,
								    2.);
	} 
      }
      for(int p = 0; p < global::L; p++){
	for(int k = 0; k < global::L; k++){
	  internal::real_L2L_S1DxS1D[1](p,k) = thetis::uniform::S1D(k,
								    thetis::uniform::point1D(p,
											     global::L,
											     1.)+1.,
								    global::L,
								    2.);
	} 
      }
    }

    template<int DIM>
    void prcmptM2ML2L(){
      internal::LpDm1 = mypow(global::L,DIM-1);
      precomputeM2M<DIM>();
      precomputeL2L<DIM>();
    }
    
    template<int DIM>
    inline void anterpolation(Vecc& bufq, vecf<DIM>& v){
      Vecc bufp(global::LL);
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::M2M_S1DxS1D[0]),
		     &bufq[0],
		     &bufp[0],
		     global::L,
		     global::L,
		     internal::LpDm1);
	}else{
	  gaia::gemm(K_(internal::M2M_S1DxS1D[1]),
		     &bufq[0],
		     &bufp[0],
		     global::L,
		     global::L,
		     internal::LpDm1);
	}
        // Swap from right to left (with modulo)
        for(int i = 0; i < global::LL; i++){
	  bufq[i] = bufp[internal::M2M_permutation[i]];
	}
      }
    }

    template<int DIM>
    inline void anterpolation(Matc& BUFQ, vecf<DIM>& v){
      Matc BUFP(global::LL,NbCol(BUFQ));
      int size = internal::LpDm1*NbCol(BUFQ);
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::M2M_S1DxS1D[0]),
		     K_(BUFQ),
		     K_(BUFP),
		     global::L,
		     global::L,
		     size);
	}else{
	  gaia::gemm(K_(internal::M2M_S1DxS1D[1]),
		     K_(BUFQ),
		     K_(BUFP),
		     global::L,
		     global::L,
		     size);
	}
        // Swap from right to left (with modulo)
        for(int k = 0; k < NbCol(BUFQ); k++){
	  cplx* BQ = K_(BUFQ)+global::LL*k;
	  cplx* BP = K_(BUFP)+global::LL*k;
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    BQ[i] = BP[internal::M2M_permutation[i]];
	  }
	}
      }
    }

    template<int DIM>
    inline void realAnterpolation(Matc& BUFQ, vecf<DIM>& v){
      Matf rBUFQ(global::LL,NbCol(BUFQ)*2);
      Matf rBUFP(global::LL,NbCol(BUFQ)*2);
      int size = internal::LpDm1*NbCol(BUFQ)*2;
      int N    = size * global::L / 2;
      // Deinterleave
      flt *prBUFQ = K_(rBUFQ);
      cplx *pBUFQ = K_(BUFQ);
#pragma omp simd
      for(int i = 0; i < N; i++){
	prBUFQ[i] = pBUFQ[i].real();
	prBUFQ[i+N] = pBUFQ[i].imag();
      }
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::real_M2M_S1DxS1D[0]),
		     K_(rBUFQ),
		     K_(rBUFP),
		     global::L,
		     global::L,
		     size);
	}else{
	  gaia::gemm(K_(internal::real_M2M_S1DxS1D[1]),
		     K_(rBUFQ),
		     K_(rBUFP),
		     global::L,
		     global::L,
		     size);
	}
        // Swap from right to left (with modulo)
        for(int k = 0; k < NbCol(rBUFQ); k++){
	  flt* BQ = K_(rBUFQ)+global::LL*k;
	  flt* BP = K_(rBUFP)+global::LL*k;
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    BQ[i] = BP[internal::M2M_permutation[i]];
	  }
	}
      }
      // Interleave
      prBUFQ = K_(rBUFQ);
      pBUFQ = K_(BUFQ);
#pragma omp simd
      for(int i = 0; i < N; i++){
        pBUFQ[i] = cplx(prBUFQ[i], prBUFQ[i+N]);
      }
    }
    

    template<int DIM>
    inline void anterpolation(cplx* bufq, vecf<DIM>& v){
      Vecc bufpp(global::LL);
      cplx* bufp = K_(bufpp);
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::M2M_S1DxS1D[0]),
		     bufq,
		     bufp,
		     global::L,
		     global::L,
		     internal::LpDm1);
	}else{
	  gaia::gemm(K_(internal::M2M_S1DxS1D[1]),
		     bufq,
		     bufp,
		     global::L,
		     global::L,
		     internal::LpDm1);
	}
        // Swap from right to left (with modulo)
        for(int i = 0; i < global::LL; i++){
	  bufq[i] = bufp[internal::M2M_permutation[i]];
	}
      }
    }
    
    template<int DIM>
    inline void interpolation(Vecc& bufq, vecf<DIM>& v){
      Vecc bufp(global::LL);
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
	if(v[DIM-1-iterIndex] < 0.){
	    gaia::gemm(K_(internal::L2L_S1DxS1D[0]),
		       &bufq[0],
		       &bufp[0],
		       global::L,
		       global::L,
		       internal::LpDm1);
	  }else{
	    gaia::gemm(K_(internal::L2L_S1DxS1D[1]),
		       &bufq[0],
		       &bufp[0],
		       global::L,
		       global::L,
		       internal::LpDm1);
	}
        // Swap from right to left (with modulo)
        for(int i = 0; i < global::LL; i++){
	  bufq[i] = bufp[internal::L2L_permutation[i]];
	}
      }
    }

    template<int DIM>
    inline void interpolation(Matc& BUFQ, vecf<DIM>& v){
      Matc BUFP(global::LL,NbCol(BUFQ));
      int size = internal::LpDm1*NbCol(BUFQ);
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::L2L_S1DxS1D[0]),
		     K_(BUFQ),
		     K_(BUFP),
		     global::L,
		     global::L,
		     size);
	}else{
	  gaia::gemm(K_(internal::L2L_S1DxS1D[1]),
		     K_(BUFQ),
		     K_(BUFP),
		     global::L,
		     global::L,
		     size);
	}
        // Swap from right to left (with modulo)
	for(int k = 0; k < NbCol(BUFQ); k++){
	  cplx* BQ = K_(BUFQ)+global::LL*k;
	  cplx* BP = K_(BUFP)+global::LL*k;
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    BQ[i] = BP[internal::L2L_permutation[i]];
	  }
	}
      }
    }

    
    template<int DIM>
    inline void realInterpolation(Matc& BUFQ, vecf<DIM>& v){
      Matf rBUFQ(global::LL,NbCol(BUFQ)*2);
      Matf rBUFP(global::LL,NbCol(BUFQ)*2);
      int size = internal::LpDm1*NbCol(BUFQ)*2;
      int N    = size * global::L / 2;
      // Deinterleave
      flt *prBUFQ = K_(rBUFQ);
      cplx *pBUFQ = K_(BUFQ);
#pragma omp simd
      for(int i = 0; i < N; i++){
	prBUFQ[i] = pBUFQ[i].real();
	prBUFQ[i+N] = pBUFQ[i].imag();
      }
      // Apply M2M
      for(int iterIndex = 0; iterIndex < DIM; iterIndex++){
	// For each 1D offset, eval 1D polynomials on 1D entries.
	// Since the same matrix is applied to a group of vectors,
	// this group is seen as a matrix => BLAS 3.
        if(v[DIM-1-iterIndex] < 0.){
	  gaia::gemm(K_(internal::real_L2L_S1DxS1D[0]),
		     K_(rBUFQ),
		     K_(rBUFP),
		     global::L,
		     global::L,
		     size);
	}else{
	  gaia::gemm(K_(internal::real_L2L_S1DxS1D[1]),
		     K_(rBUFQ),
		     K_(rBUFP),
		     global::L,
		     global::L,
		     size);
	}
        // Swap from right to left (with modulo)
        for(int k = 0; k < NbCol(rBUFQ); k++){
	  flt* BQ = K_(rBUFQ)+global::LL*k;
	  flt* BP = K_(rBUFP)+global::LL*k;
#pragma omp simd
	  for(int i = 0; i < global::LL; i++){
	    BQ[i] = BP[internal::L2L_permutation[i]];
	  }
	}
      }
      // Interleave
      prBUFQ = K_(rBUFQ);
      pBUFQ = K_(BUFQ);
#pragma omp simd
      for(int i = 0; i < N; i++){
        pBUFQ[i] = cplx(prBUFQ[i], prBUFQ[i+N]);
      }
    }

    
  } // FASTLAGRANGE
  
} // DEFMM


#endif
