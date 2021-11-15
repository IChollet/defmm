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

#ifndef TITAN_DISCRETE_FOURIER_TRANSFORMS_HPP
#define TITAN_DISCRETE_FOURIER_TRANSFORMS_HPP

#include <fftw3.h>
#include "../gaia/gaia.hpp"

namespace titan{

  namespace fft{

    template<int dim, int L>
    inline void forward(Vec<cplx>& vectorToBeTransformed, Vec<cplx>& res){
      int Larray[dim]; for(int i = 0; i < dim; i++){Larray[i] = L;}
      fftw_plan p = fftw_plan_dft(dim,
				  Larray,
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  reinterpret_cast<fftw_complex*>(K_(res)),
				  FFTW_FORWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    template<int dim, int L>
    inline void backward(Vec<cplx>& vectorToBeTransformed, Vec<cplx>& res){
      int Larray[dim]; for(int i = 0; i < dim; i++){Larray[i] = L;}
      fftw_plan p = fftw_plan_dft(dim,
				  Larray,
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  reinterpret_cast<fftw_complex*>(K_(res)),
				  FFTW_BACKWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      res /= cplx(Size(res),0.);
    }
    template<int dim, int L>
    inline void forward(Vec<cplx>& vectorToBeTransformed){
      int Larray[dim]; for(int i = 0; i < dim; i++){Larray[i] = L;}
      fftw_plan p = fftw_plan_dft(dim,
				  Larray,
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  FFTW_FORWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    template<int dim, int L>
    inline void backward(Vec<cplx>& vectorToBeTransformed){
      int Larray[dim]; for(int i = 0; i < dim; i++){Larray[i] = L;}
      fftw_plan p = fftw_plan_dft(dim,
				  Larray,
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  reinterpret_cast<fftw_complex*>(K_(vectorToBeTransformed)),
				  FFTW_BACKWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      vectorToBeTransformed /= cplx(Size(vectorToBeTransformed),0.);
    }
    inline void forward(cplx* vectorToBeTransformed, cplx* res, int L){
      fftw_plan p = fftw_plan_dft(1,
				  &L,
				  reinterpret_cast<fftw_complex*>(vectorToBeTransformed),
				  reinterpret_cast<fftw_complex*>(res),
				  FFTW_FORWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    inline void backward(cplx* vectorToBeTransformed, cplx* res, int L){
      fftw_plan p = fftw_plan_dft(1,
				  &L,
				  reinterpret_cast<fftw_complex*>(vectorToBeTransformed),
				  reinterpret_cast<fftw_complex*>(res),
				  FFTW_BACKWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      for(int i = 0; i < L; i++){
	res[i] /= cplx(double(L),0.);
      }
    }
    template<int DIM>
    inline void forward(cplx* vectorToBeTransformed, cplx* res, int L){
      int dims[DIM];
      for(int i = 0; i < DIM; i++){dims[i] = L;}
      fftw_plan p = fftw_plan_dft(DIM,
				  dims,
				  reinterpret_cast<fftw_complex*>(vectorToBeTransformed),
				  reinterpret_cast<fftw_complex*>(res),
				  FFTW_FORWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    template<int DIM>
    inline void backward(cplx* vectorToBeTransformed, cplx* res, int L){
      int dims[DIM];
      for(int i = 0; i < DIM; i++){dims[i] = L;}
      fftw_plan p = fftw_plan_dft(DIM,
				  dims,
				  reinterpret_cast<fftw_complex*>(vectorToBeTransformed),
				  reinterpret_cast<fftw_complex*>(res),
				  FFTW_BACKWARD,
				  FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      int LL = 1;
      for(int i = 1; i <= DIM; i++){LL *= L;}
      for(int i = 0; i < LL; i++){
	res[i] /= cplx(LL,0.);
      }
    }
    
    // The next functions explicitly use the column major storage of matrices.
    template<int DIM, int P>
    inline void forward(Matc& matrixToBeTransformed, Matc& res){
      int Larray[DIM]; for(int i = 0; i < DIM; i++){Larray[i] = P;}
      fftw_plan p = fftw_plan_many_dft(DIM,
				       Larray,
				       NbCol(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(res)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       FFTW_FORWARD,
				       FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    template<int DIM, int P>
    inline void backward(Matc& matrixToBeTransformed, Matc& res){
      int Larray[DIM]; for(int i = 0; i < DIM; i++){Larray[i] = P;}
      fftw_plan p = fftw_plan_many_dft(DIM,
				       Larray,
				       NbCol(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(res)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       FFTW_BACKWARD,
				       FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
      res /= cplx(NbRow(res),0.);
    }
    template<int DIM, int P>
    inline void forward(Matc& matrixToBeTransformed){
      int Larray[DIM]; for(int i = 0; i < DIM; i++){Larray[i] = P;}
      std::cout << "FORWARD :" << std::endl;
      flt t0 = gaia::chronos::time();
      fftw_plan p = fftw_plan_many_dft(DIM,
				       Larray,
				       NbCol(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       FFTW_FORWARD,
				       FFTW_ESTIMATE);
      flt t1 = gaia::chronos::time();
      fftw_execute(p);
      flt t2 = gaia::chronos::time();
      std::cout << "Plan : " << f_green << t1 - t0 << f_def << "\t" << "Execute : " << f_green << t2 - t1 << f_def << std::endl;
      fftw_destroy_plan(p);
    }
    template<int DIM, int P>
    inline void backward(Matc& matrixToBeTransformed){
      int Larray[DIM]; for(int i = 0; i < DIM; i++){Larray[i] = P;}
      std::cout << "BACKWARD :" << std::endl;
      flt t0 = gaia::chronos::time();
      fftw_plan p = fftw_plan_many_dft(DIM,
				       Larray,
				       NbCol(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       reinterpret_cast<fftw_complex*>(K_(matrixToBeTransformed)),
				       Larray,
				       1,
				       NbRow(matrixToBeTransformed),
				       FFTW_BACKWARD,
				       FFTW_ESTIMATE);
      flt t1 = gaia::chronos::time();
      fftw_execute(p);
      flt t2 = gaia::chronos::time();
      std::cout << "Plan : " << f_green << t1 - t0 << f_def << "\t" << "Execute : " << f_green << t2 - t1 << f_def << std::endl;
      fftw_destroy_plan(p);
      t0 = gaia::chronos::time();
      matrixToBeTransformed /= cplx(NbRow(matrixToBeTransformed),0.);
      t1 = gaia::chronos::time();
      std::cout << "Rescale : " << f_green << t1-t0 << f_def << std::endl;
    }
  } // fft
  
} // titan



#endif
