#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 16;
  int Ny = 5;
  double h = 4.0;

  // write results to file
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
  auto sites = SpinHalf(N);

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});

  // autompo hamiltonian
  for(auto j : lattice){
      ampo += -4, "Sz", j.s1, "Sz", j.s2;
  }
  for(auto j : range1(N)){
      ampo += -2.0*h, "Sx", j;
  }

  //initial state
  auto state = InitState(sites); 
  for(auto j : range1(N)){
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
  }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-10;

  //DMRG to find ground state at t=0
  auto H = toMPO(ampo);
  auto [energy,psi0] = dmrg(H,MPS(state),sweeps,{"Silent=",true});

  // calculate initial local energy density
  std::vector<double> LocalEnergyDM(N,0.0); // local energy density vector
  std::vector<double> LocalEnergyFit(N,0.0);
  
  //make 2D vector of ITensor for local energy operators
  std::vector<std::vector<ITensor>> LED(Nx, std::vector<ITensor>(Ny));
  std::vector<std::vector<ITensor>> LED_LR(Nx, std::vector<ITensor>(Ny));
  for(int i=1; i<=Nx; i++){ //vertical energy density
    for(int j=1; j<=Ny; j++){
      if(j==Ny){ //y-periodic boundary equations
        LED[i-1][j-1] = -4.0*sites.op("Sz",(i-1)*Ny+j)*sites.op("Sz",(i-1)*Ny+1);
      } else{
        LED[i-1][j-1] = -4.0*sites.op("Sz",(i-1)*Ny+j)*sites.op("Sz",(i-1)*Ny+j+1);
      }
    }
  }
  for(int i=1; i<Nx; i++){ //horizontal energy density
    for(int j=1; j<=Ny; j++){
      LED_LR[i-1][j-1] = -4.0*sites.op("Sz",(i-1)*Ny+j)*sites.op("Sz",i*Ny+j);
    }
  }
  
  // calculate local energy density
  for(int i=1; i<=Nx; i++){
    for(int j=1; j<=Ny; j++){ //this order to make orthogonality center easier to compute
      int index = (i-1)*Ny + j;
      psi0.position(index);
      ITensor ket;
      if(j==Ny){ //y-periodic boundary equations
        ket = psi0(index)*psi0(index-Ny+1);
        LocalEnergyDM[index-1] += elt(dag(prime(ket,"Site"))*LED[i-1][j-1]*ket);
      } else{
        ket = psi0(index)*psi0(index+1);
        LocalEnergyDM[index-1] += elt(dag(prime(ket,"Site"))*LED[i-1][j-1]*ket);
      }
      if(i<Nx){
        ket = psi0(index)*psi0(index+Ny);
        LocalEnergyDM[index-1] += elt(dag(prime(ket,"Site"))*LED_LR[i-1][j-1]*ket);
      }
      printfln("(%d,%d) = %0.3f", i, j, LocalEnergyDM[index-1]);
    }
  }
  
  // create vectors of time
  auto tval = std::vector<double>(1);
  tval.at(0) = 0.0;
  int Nt = 10;
  double dt = 0.1;
  enerfile1 << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;
  enerfile2 << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;
  enerfile1 << tval[0] << " " << energy << " " << maxLinkDim(psi0) << " "; //print to file
  enerfile2 << tval[0] << " " << energy << " " << maxLinkDim(psi0) << " "; //print to file
  for(int j = 0; j<N; j++){ //save local energy values
    enerfile1 << LocalEnergyDM[j] << " ";
    enerfile2 << LocalEnergyDM[j] << " ";
  }
  enerfile1 << std::endl;
  enerfile2 << std::endl;


  // time evolution parameters 
  auto args_DM = Args("Method=","DensityMatrix","Cutoff=",1E-10,"MaxDim=",3000);
  auto args_Fit = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",3000);

  //time evolution fields h and diff
  auto hval = std::vector<double>(1);
  auto diff = std::vector<double>(1);
  double step = (h-2.0)/Nt;
  hval.at(0) = h-step;
  diff.at(0) = -step;
  for(int j=2; j<=Nt; j++){ //make linear ramp
    hval.push_back(h-step*j);
    diff.push_back(-step);
  }
  
  // initial conditions
  auto psi_DM = psi0; //keep psi0 for future reference
  auto psi_Fit = psi0;
  
  //
  // time evolve
  //
  for(int i=0; i<Nt; i++){
    //
    // update autoMPO
    for(auto j : range1(N)){
      ampo += -2.0*diff[i], "Sx", j;
    }
    //
    //time evolution operators
    auto expH1 = toExpH(ampo, 0.5*dt*(1+Cplx_i)); //time evolve by 0.5*(1+i)*dt
    auto expH2 = toExpH(ampo, 0.5*dt*(1-Cplx_i)); //time evolve by 0.5*(i-i)*dt

    tval.push_back(tval[i]+dt); //update time vector
    // check DensityMatrix method for MPO*MPS
    psi_DM = applyMPO(expH1,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_DM = applyMPO(expH2,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels

    // calculate local energy density
    for(auto b : lattice){
      psi_DM.position(b.s1);
      auto ket = psi_DM(b.s1)*psi_DM(b.s2);
      auto LED = -4.0*sites.op("Sz",b.s1)*sites.op("Sz",b.s2);
      LocalEnergyDM[b.s1] += inner(ket,LED,ket);
    }
    for(int j = 0; j<N; j++){ //save local energy values
    enerfile1 << LocalEnergyDM[j] << " ";
    }
    enerfile1 << std::endl;

    //check Fit method for MPO*MPS
    psi_Fit = applyMPO(expH1,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_Fit = applyMPO(expH2,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels

    // calculate local energy density
    for(auto b : lattice){
      psi_DM.position(b.s1);
      auto ket = psi_Fit(b.s1)*psi_Fit(b.s2);
      auto LED = -4.0*sites.op("Sz",b.s1)*sites.op("Sz",b.s2);
      LocalEnergyFit[b.s1] += inner(ket,LED,ket);
    }
    for(int j = 0; j<N; j++){ //save local energy values
    enerfile2 << LocalEnergyFit[j] << " ";
    }
    enerfile2 << std::endl;
    
    auto energy_DM = innerC(psi_DM,H,psi_DM).real();
    auto energy_Fit = innerC(psi_Fit,H,psi_Fit).real();

    enerfile1 << tval[i+1] << " " << energy_DM << " " << maxLinkDim(psi_DM) << std::endl;
    enerfile2 << tval[i+1] << " " << energy_Fit << " " << maxLinkDim(psi_Fit) << std::endl;
    printfln("Iteration %d, DensityMatrix energy = %0.3g, max link dim is %d",i+1,energy_DM,maxLinkDim(psi_DM));
    printfln("            ,          Fit energy = %0.3g, max link dim is %d",energy_Fit,maxLinkDim(psi_Fit));

  }
  
  enerfile1.close(); enerfile2.close();

  return 0;
  }
