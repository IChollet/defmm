#include "../include/interface.hpp"

using namespace defmm;

int main(int argc, char* argv[]){

  const int  DIM   = 3;
  const int  ORDER = (argc > 1 ? atoi(argv[1]) : 6);
  const int  NCRIT = (argc > 2 ? atoi(argv[2]) : 32);
  const flt  KAPPA = 1.;
  
  // FMM matrix
  IBFMM_Mat<DIM> A;

  // Add point clouds
  int N = 41334;
  A.addSourceParticlesINP("./mesh/refinedcube.inp",N);
  A.addTargetParticlesINP("./mesh/refinedcube.inp",N);
  std::cout << "Refined cube test case with N="  << f_cyan << N << f_def << std::endl;
  std::cout << "kappa D=64" << std::endl;

  // Get random charge vector
  Vecc In(N), Out(N); Out = cplx0;
  for(int n = 0; n < N; n++){
    In[n] = cplx(urand);}

  // Precompute
  A.prcmpt(ORDER,NCRIT,KAPPA);
  
  // Apply FMM
  gemv(A,In,Out);

  // Comparison with direct computation
  CheckFULL(A);
  
  
  return 0;
}
