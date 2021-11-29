#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 64;
  int Ny = 5;
  double h=4.0;

  //write results to file
  char schar1[50], schar2[50];
  int n1 = std::sprintf(schar1,"Nx_%d_Ny_%d_h_%1.3g_Ising2d_teval_MPO-DM.dat",Nx,Ny,h);
  int n2 = std::sprintf(schar2,"Nx_%d_Ny_%d_h_%1.3g_Ising2d_teval_MPO-Fit.dat",Nx,Ny,h);
  std::string s1(schar1);
  std::string s2(schar2);
  std::ofstream enerfile1;
  std::ofstream enerfile2;
  enerfile1.open(s1); // opens the file
   if( !enerfile1 ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
   }
  enerfile2.open(s2); // opens the file
   if( !enerfile2 ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
   }

  auto N = Nx * Ny;
  auto sites = SpinHalf(N,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});

  // autompo hamiltonian
  for(auto j : lattice){
      ampo += -4, "Sz", j.s1, "Sz", j.s2;
      }
  for(auto j : range1(N)){
      ampo += -2.0*h, "Sx", j;
      }
  
  auto state = InitState(sites); //initial state
  for(auto j : range1(N)){
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
      }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-8;

  //DMRG to find ground state at t=0
  auto H = toMPO(ampo);
  auto [energy,psi0] = dmrg(H,MPS(state),sweeps,{"Silent=",true});
  //PrintData(H);
  
  // create vectors of time
  auto tval = std::vector<double>(1);
  tval.at(0) = 0.0;
  int Nt = 10;
  double dt = 0.1;
  enerfile1 << "tval" << " " << "energy" << " " << "MaxDim" << std::endl;
  enerfile1 << tval[0] << " " << energy << " " << maxLinkDim(psi0) << std::endl; //print to file

  //time evolution operators
  auto expH1 = toExpH(ampo, 0.5*dt*(1+Cplx_i)); //time evolve by 0.5*(1+i)*dt
  auto expH2 = toExpH(ampo, 0.5*dt*(1-Cplx_i)); //time evolve by 0.5*(i-i)*dt
  auto args_DM = Args("Method=","DensityMatrix","Cutoff=",1E-8,"MaxDim=",3000);
  auto args_Fit = Args("Method=","Fit","Cutoff=",1E-8,"MaxDim=",3000);
  //PrintData(expH1);
  //PrintData(expH2);

  //time evolution fields h and diff
  auto hval = std::vector<double>(1);
  auto diff = std::vector<double>(1);
  double step = (h-2.0)/Nt;
  for(int j=1; j<=Nt; j++){ //make linear ramp
    hval.push_back(h-step*j);
    diff.push_back(-step);
    println(hval[j-1]);
    println(diff[j-1]);
  }
  
  /*
  auto psi_DM = psi0; //keep psi0 for future reference
  auto psi_Fit = psi0;
  // time evolve
  for(int i=0; i<Nt; i++){
    for(auto j : range1(N)){
      ampo += -2.0*diff[i], "Sx", j;
      }
    tval.push_back(tval[i]+dt);
    // check DensityMatrix method for MPO*MPS
    psi_DM = applyMPO(expH1,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_DM = applyMPO(expH2,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels

    //check Fit method for MPO*MPS
    psi_Fit = applyMPO(expH1,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_Fit = applyMPO(expH2,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels
    
    auto energy_DM = innerC(psi_DM,H,psi_DM).real();
    auto energy_Fit = innerC(psi_Fit,H,psi_Fit).real();

    enerfile1 << tval[i+1] << " " << energy_DM << " " << maxLinkDim(psi_DM) << std::endl;
    enerfile2 << tval[i+1] << " " << energy_Fit << " " << maxLinkDim(psi_Fit) << std::endl;
    printfln("Iteration %d, DensityMatrix energy = %0.3g, max link dim is %d",i+1,energy_DM,maxLinkDim(psi_DM));
    printfln("            ,          Fit energy = %0.3g, max link dim is %d",energy_Fit,maxLinkDim(psi_Fit));

  }
  */
  enerfile1.close();

  return 0;
  }
