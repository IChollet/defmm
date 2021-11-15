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

#ifndef THETIS_KERNELS_HEADER_HPP
#define THETIS_KERNELS_HEADER_HPP

#include "thetis_types.hpp"

namespace thetis{
  typedef cplx(*kernel_op)(const flt&,const flt&);
  namespace kernels{
    inline cplx Helmholtz(const flt& R, const flt& kappa){return exp(cplxi*kappa*R)/(pi4*R);}
    inline cplx Laplace  (const flt& R, const flt& kappa){return 1./(pi4*R);}
    inline cplx Gaussian (const flt& R, const flt& kappa){return exp(-R*R);}
    inline cplx MultiFreq(const flt& R, const flt& kappa){return exp(cplxi*kappa*R)/(pi4*R)+exp(cplxi*kappa*R/3.6)/(pi4*R)+exp(cplxi*kappa*R/2.)/(pi4*R)+exp(cplxi*kappa*R/8.)/(pi4*R);}
    inline cplx Identity (const flt& R, const flt& kappa){return R;}
  } // KERNELS
}

#endif
