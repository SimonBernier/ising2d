#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 32;
  int Ny = 6;
  double h=4.0;

  //write results to file
  char schar1[50];
  int n1 = std::sprintf(schar1,"Nx_%d_Ny_%d_h_%1.3g_Ising2d_teval.dat",Nx,Ny,h);
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

  // autompo hamiltonian
  for(auto j : lattice)
      {
      ampo += -1, "Sz", j.s1, "Sz", j.s2;
      }
  for(auto j : range1(N))
      {
      ampo += -h, "Sx", j;
      }
  
  auto state = InitState(sites); //initial state
  for(auto j : range1(N))
      {
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-8;

  //DMRG to find ground state at t=0
  auto H = toMPO(ampo);
  auto [energy,psi] = dmrg(H,state,sweeps,{"Silent=",true});

  // create vectors of time
  auto tval = std::vector<double>(1);
  tval.at(0) = 0.0;
  int Nt = 10;
  double dt = 0.1;
  enerfile1 << "tval" << " " << "energy" << std::endl;
  enerfile1 << tval[0] << " " << energy << std::endl; //print to file

  //time evolution operators
  auto expH1 = toExpH(ampo, 0.5*dt*(1+Cplx_i)); //time evolve by 0.5*(1+i)*dt
  auto expH2 = toExpH(ampo, 0.5*dt*(1-Cplx_i)); //time evolve by 0.5*(i-i)*dt
  auto args = Args("Method=","DensityMatrix","Cutoff=",1E-9,"MaxDim=",3000);

  // time evolve
  for(int i=0; i<Nt; i++){
    tval.push_back(tval[i]+dt);
    psi = applyMPO(expH1,psi,args);
    psi = applyMPO(expH2,psi,args);
    psi.noPrime().normalize(); //need to do this after each?
    
    energy = inner(psi,H,psi);

    enerfile1 << tval[i+1] << " " << energy << std::endl;
    printfln("Iteration %d, energy = %0.3g",i,energy);

  }
  enerfile1.close();

  return 0;
  }
