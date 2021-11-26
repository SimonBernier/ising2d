#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 64;
  int Ny = 10;

  //write results to file
  char schar1[50];
  int n1 = std::sprintf(schar1,"Nx_%d_Ny_%d_Ising2d_Gap.dat",Nx,Ny);
  std::string s1(schar1);
  std::ofstream enerfile1;
  enerfile1.open(s1); // opens the file
   if( !enerfile1 ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
   }

  auto N = Nx * Ny;
  auto sites = SpinHalf(N,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});

  // create vectors of h
  auto h = std::vector<double>(1);
  h.at(0) = 1.0;
  int iter = 10;
  double h_step = 0.05;
  for(int i=1; i<=iter; i++){
      h.push_back(h[0] + i*h_step);
  }
  auto diff = std::vector<double>(1); //used to create new Hamiltonian
  diff.at(0) = 0.0;
  for(int i=1; i<=iter; i++){
      diff.push_back(h_step);
  }
  // autompo hamiltonian
  for(auto j : lattice)
      {
      ampo += -4, "Sz", j.s1, "Sz", j.s2;
      }
  for(auto j : range1(N))
      {
      ampo += -2.0*h[0], "Sx", j;
      }
  
  auto state = InitState(sites); //initial state
  for(auto j : range1(N))
      {
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  // 2d ising model parameters
  auto sweeps = Sweeps(8);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-8;

  enerfile1 << "hval" << " " << "maxBondDimGS" << " " << "maxBondDimExc" << " " << "energy" << " " << "gap" << " " << std::endl;
  
  // loop over values of h
  for(int i=0; i<iter; i++){
    for(auto j : range1(N)){
        ampo += 2*diff[i], "Sx", j;
        }
    auto H = toMPO(ampo); //12x12 matrices
    auto [energy,psi0] = dmrg(H,state,sweeps,{"Silent=",true});

    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;

    //
    // Here the Weight option sets the energy penalty for
    // psi1 having any overlap with psi0
    //
    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    auto gap = en1-energy; //compute gap energy

    enerfile1 << h[i] << " " << maxLinkDim(psi0) << " " << maxLinkDim(psi1) << " " << energy << " " << gap << " " << std::endl;
    printfln("Iteration %d, h = %.3g, gap = %0.3g",i,h[i],gap);

  }
  enerfile1.close();

  return 0;
  }
