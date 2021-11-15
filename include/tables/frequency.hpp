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

#ifndef FREQUENCY_MAIN_HEADER_HPP
#define FREQUENCY_MAIN_HEADER_HPP

#include "../structs/tree.hpp"
#include "../global.hpp"
#include "../symmetries/permutations.hpp"
#include "newdirections.hpp"

namespace defmm{

  /*
    Decides if a stage is lf or hf.
    If a stage is lf, the the corresponding entry
    in boolArrayHFStage is 0.
    If it is HF, the corresponding entry is its 
    level in the direction tree.
  */
  template<int DIM>
  inline void fixFrequencyRegime(Tree<DIM>& T,
				 Tree<DIM>& S,
				 flt lowFrequencyThreshold){
    int nlvl = std::max(Nlvl(T),Nlvl(S));
    // Find first HF level
    flt stageDiameter = 2.*Radius(T);
    int cpt           = 0;
    if(stageDiameter*global::kappa > lowFrequencyThreshold){
      while(stageDiameter*global::kappa > lowFrequencyThreshold){
        stageDiameter /= 2.;
	cpt++;
      }
    }else{
      cpt = 0;
    }
    // Get the associated direction tree
    if(cpt){iterativeGetDirections(global::directions, cpt);}
    // Test if there exist a HF regime
    if(cpt > 1){
      // There exist a high frequency regime
      global::noHF = false;
    }else{
      // There is no hf regime
      global::noHF = true ;
    }
    // Set HF stages to their HF level
    int hfLvl = 1;
    for(int l = cpt; l > -1; l--){
      if(l < nlvl){
	global::isHFStage[l] = hfLvl;
	global::HFStage[l] = hfLvl-1;
      }
      hfLvl++;
    }
    for(int l = nlvl-1; l > cpt; l--){
      global::isHFStage[l] = 0;
      global::HFStage[l]   = 0;
    }
  }

  /*
    Fix the number of direction (modulo permutations)
    per hf level, using the convention 0 for low frequency
    stages.
  */
  template<int DIM>
  void fixNumberOfDirectionPerLevel(Tree<DIM>& T, Tree<DIM>& S){
    for(int l = std::max(Nlvl(T),Nlvl(S))-1; l > -1; l--){
      if(global::isHFStage[l] && l != 0){
        global::nDirPerStage[l] = global::directions.indices[global::HFStage[l]];
      }else{
	global::nDirPerStage[l] = 1;
      }
    }
  }

  /*
    Upper stages in the high frequency regime are empty
    in term of interaction (because the distances between
    two cells of those stages are not admissible w.r.t. 
    the directional admissibility condition).
    This function associates to each level a boolean which
    is true iff their is a possible interaction on the stage.
  */
  template<int DIM>
  void fixPrecomputationPerLevel(Tree<DIM>& T, Tree<DIM>& S){
    int maxStage = 0;
    flt BoxDiameter = 2.*Radius(T);
    while( 1. < global::kappa*BoxDiameter/flt(1 << maxStage)/flt(1 << maxStage) ){
      maxStage++;
      if(maxStage >= std::max(Nlvl(T),Nlvl(S))){
	std::cout << std::max(Nlvl(T),Nlvl(S)) << "\t" << maxStage << std::endl; 
	std::cout << f_red << "Error in 'frequency.hpp'" << std::endl
		  <<          "   -> In 'fixPrecomputationPerLevel'" << std::endl
		  <<          "      -> No precomputation allowed!" << f_def << std::endl;
	exit(1);
      }
    }
    maxStage = std::max(maxStage,2);
    for(int l = maxStage; l < std::max(Nlvl(T),Nlvl(S)); l++){
      global::precomputationAtStage[l] = true;
    }
    for(int l = 0; l < maxStage; l++){
      global::precomputationAtStage[l] = false;
    }
  }

  /*
    Fix lower bounds for M2L tables.
    All bounds are the same in the low frequency regime.
    In the HF regime, the pattern depends on the level.
  */
  template<int DIM>
  void fixTablesLowerBounds(Tree<DIM>& T,
			    Tree<DIM>& S,
			    int*       boolArrayHFStage,
			    bool*      precomputeAtStage,
			    int*       tableLowerBounds){
    // All fixed at '2' (as in the lf case).
    // The HF precomputation will automatically
    // detect the effective M2L that need to
    // be computed, using HF-admissibility.
    // This needs to be ehanced later...
    for(int l = 0; l < std::max(Nlvl(T),Nlvl(S)); l++){
      tableLowerBounds[l] = 2;
    }
  }

  /*
    Fix upper bounds for M2L tables.
    Once again, those bounds depend on the frequency
    regime. They have to be balanced by the maximum 
    number of leaves per stage.
  */
  template<int DIM>
  void fixTablesUpperBounds(Tree<DIM>& T,
			    Tree<DIM>& S,
			    flt        eta,
			    int*       boolArrayHFStage,
			    bool*      precomputeAtStage,
			    int*       tableUpperBounds){
    for(int l = 0; l < std::max(Nlvl(T),Nlvl(S)); l++){
      if(precomputeAtStage[l]){
	if(boolArrayHFStage[l]){
	  // HF case
	  int nbCellOnStage   = (1<<l);
	  tableUpperBounds[l] = nbCellOnStage;
	}else{
	  // LF case
	  // (Next test allowed because precomputationAtStage[0] = false)
	  if(global::isHFStage[l-1]){
	    tableUpperBounds[l] = 4*std::max(2,int(eta));
	  }else{
	    tableUpperBounds[l] = 4;
	  }
	}
      }else{
	tableUpperBounds[l] = -1;
      }
    }
  }

  /*
    Print infos about frequency regimes on different levels.
  */
  template<int DIM>
  void printFrequencyStageInfos(Tree<DIM>& T, Tree<DIM>& S){
    std::cout << "FREQUENCY infos" << std::endl;
    std::cout << "\t[HF level, Nb directions, Upper bound, Min dist]" << std::endl;
    for(int l = 0; l < std::max(Nlvl(T),Nlvl(S)); l++){
      std::cout << "\tLevel : "      << f_cyan << l << f_def << std::endl
		<< "\t\tNb lambda     : " << f_cyan << 2.*Radius(T)/(1 << l)*global::kappa << f_def << std::endl;
      if(global::precomputationAtStage[l]){
	std::cout << "\t\t[" << f_cyan << global::isHFStage[l] << f_def
		  << ", " << f_cyan << global::nDirPerStage[l] << f_def
		  << ", " << f_cyan << global::nCellPerDirMax[l] << f_def
		  << ", " << f_cyan << std::max(pow(2.*Radius(T)/(1 << l),2)*global::kappa, 2.*Radius(T)/(1 << l)) << f_def << "]" << std::endl;
      }
    }
  }

  /*
    Allocation for all arrays of frequency.hpp
  */
  template<int DIM>
  void allocateFrequency(Tree<DIM>& T, Tree<DIM>& S){
    int nlvl = std::max(Nlvl(T),Nlvl(S));
    global::isHFStage             = new int [nlvl];
    global::precomputationAtStage = new bool[nlvl];
    global::nCellPerDirMin        = new int [nlvl];
    global::nCellPerDirMax        = new int [nlvl];
    global::nDirPerStage          = new int [nlvl];
    global::HFStage               = new int[nlvl];
  }

  template<int DIM>
  void initFrequency(Tree<DIM>& T, Tree<DIM>& S){
    flt eta = 1.;
    allocateFrequency(T,S);
    fixFrequencyRegime(T,S,eta);
    fixNumberOfDirectionPerLevel(T,S);
    fixPrecomputationPerLevel(T,S);
    fixTablesLowerBounds(T,
			 S,
			 global::isHFStage,
			 global::precomputationAtStage,
			 global::nCellPerDirMin);
    fixTablesUpperBounds(T,
			 S,
			 eta,
			 global::isHFStage,
			 global::precomputationAtStage,
			 global::nCellPerDirMax);
#ifdef DEFMM_VERBOSE
    printFrequencyStageInfos(T,S);
#endif
  }

} // DEFMM

#endif 
