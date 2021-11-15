#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 64;
  int Ny = 3;

  //write results to file
  char schar1[50];
  int n1 = std::sprintf(schar1,"Nx_%d_Ny_%d_Ising2d_Gap.txt",Nx,Ny);
  std::string s1(schar1);
  std::ofstream enerfile1(s1);
  if (!enerfile1.is_open())
    {
        return 0;
    }

  double h = 4.0; //initialize h in the gapped phase

  auto N = Nx * Ny;
  auto sites = SpinHalf(N,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});

  for(auto j : lattice)
      {
      ampo += -4, "Sz", j.s1, "Sz", j.s2;
      }
  for(auto j : range1(N))
      {
      ampo += -2*h, "Sx", j;
      }
  auto H = toMPO(ampo); //12x12 matrices
  
  auto state = InitState(sites);
  for(auto j : range1(N))
      {
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-8;

  std::cout << "hval" << " " << "maxBondDimGS" << " " << "maxBondDimExc" << " " << "energy" << " " << "gap" << " " << std::endl;

  auto [energy,psi0] = dmrg(H,state,sweeps,{"Quiet=",true});

  auto wfs = std::vector<MPS>(1);
  wfs.at(0) = psi0;

  //
  // Here the Weight option sets the energy penalty for
  // psi1 having any overlap with psi0
  //
  auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Quiet=",true,"Weight=",20.0});
  auto gap = en1-energy; //compute gap energy

  std::cout << h << " " << maxLinkDim(psi0) << " " << maxLinkDim(psi1) << " " << energy << " " << gap << " " << std::endl;
  enerfile1.close();

  return 0;
  }
