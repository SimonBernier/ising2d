#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 64;
  int Ny = 10;
  double h = 4.0;
  if(argc > 3)
    h = std::stof(argv[3]);
  if(argc > 2)
    Ny = std::stoi(argv[2]);
  if(argc > 1)
    Nx = std::stoi(argv[1]);

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
  //auto H = toMPO(ampo,{"Exact=",true}); //14x14 matrices

  PrintData(H);

  /*
  auto state = InitState(sites);
  for(auto j : range1(N))
      {
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  // 2d hubbard model parameters
  //auto sweeps = Sweeps(15);
  //sweeps.maxdim() = 20, 60, 100, 100, 200, 400, 800, 2000, 3000;
  //sweeps.noise() = 1E-7, 1E-8, 1E-10, 0;
  //sweeps.cutoff() = 1E-6;

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200;
  sweeps.cutoff() = 1E-8;

  PrintData(sweeps);

  auto psi0 = randomMPS(state);
  auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet=",true});

  PrintData(Nx);
  PrintData(Ny);
  PrintData(h);
  //PrintData(totalQN(psi));
  PrintData(maxLinkDim(psi));
  PrintData(energy);
  */

  return 0;
  }
