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

#ifndef THETIS_POLYNOMIALS_HPP
#define THETIS_POLYNOMIALS_HPP

#include "../gaia/gaia.hpp"

namespace thetis{

  namespace scaling{
    inline flt Phi(flt x, flt a, flt b){
      return 0.5*(a+b)+0.5*(b-a)*x;
    }
    template <int DIM>
    inline vecf<DIM> Phi(const vecf<DIM>& X, const vecf<DIM>& A, const vecf<DIM>& B){
      vecf<DIM> res; for(int d=0; d<DIM; d++){ res[d] = Phi(X[d],A[d],B[d]);}
      return res;
    }
    inline flt iPhi(const flt& x, const flt& a, const flt& b){return (2.*x-b-a)/(b-a);}
    template<int DIM>
    inline vecf<DIM> iPhi(const vecf<DIM>& X, const vecf<DIM>& A, const vecf<DIM>& B){
      vecf<DIM> res; for(int d=0; d<DIM; d++){ res[d] = iPhi(X[d],A[d],B[d]);}
      return res;
    }
  }
  
  namespace uniform{
    template<int P>
    flt S1D(int i, const flt& x, const flt& a, const flt& b){
      flt res = 1.;
      flt nx  = scaling::iPhi(x,a,b);
      flt pi  = -1.+2.*((flt)( i)/(flt)(P-1));
      for(int ii=0; ii<i; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (nx-pii)/(pi-pii);
      }
      for(int ii=i+1; ii<P; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (nx-pii)/(pi-pii);
      }
      return res;
    }
    template<int P>
    flt S1D(int i, const flt& x){
      flt res = 1.;
      flt pi  = -1.+2.*((flt)( i)/(flt)(P-1));
      for(int ii=0; ii<i; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (x-pii)/(pi-pii);
      }
      for(int ii=i+1; ii<P; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (x-pii)/(pi-pii);
      }
      return res;
    }
    
    flt S1D(int i, const flt& x, int P, const flt& alpha){
      flt res = 1.;
      flt nx  = x/alpha;
      flt pi  = -1.+2.*((flt)( i)/(flt)(P-1));
      for(int ii=0; ii<i; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (nx-pii)/(pi-pii);
      }
      for(int ii=i+1; ii<P; ii++){
	flt pii = -1.+2.*((flt)(ii)/(flt)(P-1));
	res *= (nx-pii)/(pi-pii);
      }
      return res;
    }
    
    template<int DIM, int L>
    inline flt S(const veci<DIM>& I, const vecf<DIM>& X, const vecf<DIM>& A, const vecf<DIM>& B){
      flt res = 1.;
      for(int d = 0; d < DIM; d++){
	res *= S1D<L>(I[d],X[d],A[d],B[d]);
      }
      return res;
    }
    template<int DIM, int L>
    inline flt S(const veci<DIM>& I, const vecf<DIM>& X){
      flt res = 1.;
      for(int d = 0; d < DIM; d++){
	res *= S1D<L>(I[d],X[d]);
      }
      return res;
    }
    template<int DIM>
    inline flt S(const veci<DIM>& I, const vecf<DIM>& X, int L, const flt& alpha){
      flt res = 1.;
      for(int d = 0; d < DIM; d++){
        res *= S1D(I[d],X[d],L,alpha);
      }
      return res;
    }
    template<int L> inline flt point1D(int i, const flt& a, const flt& b){
      return scaling::Phi(-1.+2.*((flt)( i)/(flt)(L-1)), a, b);}
    template<int L> inline flt point1D(int i){
      return -1.+2.*((flt)( i)/(flt)(L-1));}
    inline flt point1D(int i, int L){
      return -1.+2.*((flt)( i)/(flt)(L-1));}
    inline flt point1D(int i,int L, const flt& alpha){
      return alpha*(-1.+2.*((flt)( i)/(flt)(L-1)));}
    template<int DIM, int L>
    inline vecf<DIM> point(const veci<DIM>& I, const vecf<DIM>& A, const vecf<DIM>& B){
      vecf<DIM> res; for(int d = 0; d < DIM; d++){res[d] = point1D<L>(I[d],A[d],B[d]);}
      return res;
    }
    template<int DIM, int L>
    inline vecf<DIM> point(const veci<DIM>& I){
      vecf<DIM> res; for(int d = 0; d < DIM; d++){res[d] = point1D<L>(I[d]);}
      return res;
    }
    template<int DIM>
    inline vecf<DIM> point(const veci<DIM>& I, int L, const flt& alpha){
      vecf<DIM> res; for(int d = 0; d < DIM; d++){res[d] = point1D(I[d],L,alpha);}
      return res;
    }
  } // uniform
  
  namespace chebyshev{
    // Ce serait peut-être bien de tabuler tout ça... (comme dans dfmm)
    template<int L> inline flt point1D(int i, flt a, flt b){return scaling::Phi(cos((2.*flt(i)+1.)/flt(L)*M_PI*0.5),a,b);}
    template<int L> inline flt point1D(int i){return cos((2.*flt(i)+1.)/flt(L)*M_PI*0.5);}
    inline flt T(int n, const flt& x){return cos(flt(n)*acos(x));}
    template<int P>
    inline flt S1D(int i, const flt& x, const flt& a, const flt& b){
      flt one_over_P_plus_one = 1./flt(P); // Because we block at P-1 with '< P', not at P 
      flt res = 0.;
      flt xab = scaling::iPhi(x,a,b);
      for(int n=1; n < P; n++){
	res += T(n,xab)*T(n,point1D<P>(i));
      }
      res *= 2.*one_over_P_plus_one;
      res += one_over_P_plus_one;
      return res;
    }
    template<int DIM, int L>
    inline vecf<DIM> point(const veci<DIM>& I, const vecf<DIM>& A, const vecf<DIM>& B){
      vecf<DIM> res;
      for(int d = 0; d < DIM; d++){
	res[d] = point1D<L>(I[d],A[d],B[d]);
      }
      return res;
    }
    template<int DIM, int L>
    inline flt S(const veci<DIM>& I, const vecf<DIM>& X, const vecf<DIM>& A, const vecf<DIM>& B){
      flt res = 1.;
      for(int d = 0; d < DIM; d++){
	res *= S1D<L>(I[d],X[d],A[d],B[d]);
      }
      return res;
    }
  } // chebyshev  
  
} // thetis


#endif
