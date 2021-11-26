#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int N = 16;

  //write results to file
  char schar1[50];
  int n1 = std::sprintf(schar1,"N_%d_IsingGapScaling_10sweeps.dat",N);
  std::string s1(schar1);
  std::ofstream enerfile1;
  enerfile1.open(s1); // opens the file
   if( !enerfile1 ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
   }

  auto sites = SpinHalf(N,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);

  // create vectors of h
  auto h = std::vector<double>(1);
  h.at(0) = 2.0;
  int iter = 20;
  double h_step = 2.0/iter;
  for(int i=1; i<=iter; i++){
      h.push_back(h[0] - i*h_step);
  }
  auto diff = std::vector<double>(1); //used to create new Hamiltonian
  diff.at(0) = 0.0;
  for(int i=0; i<iter; i++){
      diff.push_back(h[i+1]-h[i]);
  }
  
  // autompo hamiltonian
  for(int j=1; j<N; j++){
      ampo += -4.0, "Sz", j, "Sz", j+1;
      }
  for(int j=1; j<=N; j++){
      ampo += -2.0*h[0], "Sx", j;
      }

  // 2d ising model parameters
  auto sweeps = Sweeps(10);
  sweeps.maxdim() = 20, 50, 100, 200;
  sweeps.cutoff() = 1E-10;

  enerfile1 << "hval" << " " << "gsEn" << " " << "excEn" << " " << "overlap" << " " << "var0" << " " << "var1" << " " << std::endl;
  
  // loop over values of h
  for(int i=0; i<=iter; i++){
    for(auto j : range1(N)){
        ampo += -2.0*diff[i], "Sx", j;
    }
    auto H = toMPO(ampo);
    auto [en0,psi0] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
    auto var0 = inner(H,psi0,H,psi0) - en0*en0;

    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;

    //
    // Here the Weight option sets the energy penalty for
    // psi1 having any overlap with psi0
    //
    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    auto var1 = inner(H,psi1,H,psi1) - en1*en1;
    auto gap = en1-en0; //compute gap energy

    enerfile1 << h[i] << " " << en0 << " " << en1 << " " << inner(psi1,psi0) << " " << var0 << " " << var1 << " " << std::endl;
    printfln("Iteration %d, h = %.3g, gap = %0.3g",i,h[i],gap);
    printfln("              var0 = %0.3e, var1 = %0.3e, overlap = %0.3e", var0, var1, inner(psi1,psi0));

    }
  enerfile1.close();
  

  return 0;
  }
