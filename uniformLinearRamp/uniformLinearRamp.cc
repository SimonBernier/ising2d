#include "itensor/all.h"

using namespace itensor;

//function definition for calculation of local energy
std::vector<double> localEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>, std::vector<ITensor>);

int main(int argc, char *argv[])
  {
  int Nx = 16;
  int Ny = 4;
  float h = 4.0;
  int Nt = 20;
  float dt = 0.1;
  float tau = 1.0;
  if(argc > 6)
    tau = std::stof(argv[6]);
  if(argc > 5)
    dt = std::stof(argv[5]);
  if(argc > 4)
    Nt = std::stoi(argv[4]);
  if(argc > 3)
    h = std::stof(argv[3]);
  if(argc > 2)
    Ny = std::stoi(argv[2]);
  if(argc > 1)
    Nx = std::stoi(argv[1]);

  // write results to file
  char schar[50];
  int n1 = std::sprintf(schar,"Nx_%d_Ny_%d_h_%0.3g_tau_%0.3g_uniLin.dat",Nx,Ny,h,tau);
  std::string s1(schar);
  std::ofstream enerfile;
  enerfile.open(s1); // opens the file
  if( !enerfile ) { // file couldn't be opened
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

  //initial state
  auto state = InitState(sites); 
  for(auto j : range1(N)){
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
  }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-10;

  // calculate initial local energy density
  std::vector<double> LocalEnergy(N,0.0); // local energy density vector
  
  //make 2D vector of ITensor for local energy operators
  //long-range interactions have the same structure as nearest-neighbour when we use swap gates
  std::vector<std::vector<ITensor>> LED(Nx, std::vector<ITensor>(Ny));
  std::vector<ITensor> LEDyPBC(Nx);
  
  for(int i=1; i<=Nx; i++){
    for(int j=1; j<=Ny; j++){
      int index = (i-1)*Ny+j;
      //MPS nearest-neighbour
      //good for long-range interaction where index+Ny is brought to index+1
      if(index<N){
        LED[i-1][j-1] = -4.0*sites.op("Sz",index)*sites.op("Sz",index+1);
      }
      //y-periodic boundary equations
      if(j==Ny){
        // site index-Ny+1 is moved to site index-1 with swap gates
        LEDyPBC[i-1] = -4.0*sites.op("Sz",index-1)*sites.op("Sz",index);
      }
    }
  }
  
  //DMRG to find ground state at t=0
  auto H = toMPO(ampo);
  auto [energy,psi] = dmrg(H,MPS(state),sweeps,{"Silent=",true});
  
  // calculate local energy density
  LocalEnergy = localEnergy(Nx, Ny, sites, psi, LED, LEDyPBC);

  //store to file
  enerfile << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;
  enerfile << 0.0 << " " << energy << " " << maxLinkDim(psi) << " "; //print to file
  for(int j = 0; j<N; j++){ //save local energy values
    enerfile << LocalEnergy[j] << " ";
  }
  enerfile << std::endl;

  // time evolution parameters
  auto args_DM = Args("Method=","DensityMatrix","Cutoff=",1E-10,"MaxDim=",3000);
  auto args_Fit = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",3000);

  //time evolution fields h and diff
  std::vector<double> tval(Nt, 0.0); //time vector
  auto hval = std::vector<double>(Nt+1,0.0);
  auto diff = std::vector<double>(Nt+1,0.0);
  int Nstep = tau/dt;
  double step = (h-2.0)/Nstep;
  hval[0] = h;
  diff[0] = 0.0;
  for(int j=1; j<=Nt; j++){ //make linear ramp
    if(j<=Nstep){
      hval[j] = h-step*j;
    } else{
      hval[j] = hval[j-1];
    }
    diff[j] = hval[j]-hval[j-1];
  }
  
  //
  // time evolve
  //
  for(int n=0; n<Nt; n++){
    //
    // update autoMPO
    for(auto j : range1(N)){
      ampo += -2.0*diff[n], "Sx", j;
    }
    //
    //time evolution operators
    auto expH1 = toExpH(ampo, 0.5*dt*(1+Cplx_i)); //time evolve by 0.5*(1+i)*dt
    auto expH2 = toExpH(ampo, 0.5*dt*(1-Cplx_i)); //time evolve by 0.5*(i-i)*dt
    tval[n]=(n+1)*dt; //update time vector

    //Fit method for MPO*MPS
    psi = applyMPO(expH1,psi,args_Fit);
    psi.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi = applyMPO(expH2,psi,args_Fit);
    psi.noPrime().normalize(); //need to do this after each to take care of prime levels
    auto energy = innerC(psi,H,psi).real();

    // calculate local energy density
    LocalEnergy = localEnergy(Nx, Ny, sites, psi, LED, LEDyPBC);

    //write to file
    enerfile << tval[n] << " " << energy << " " << maxLinkDim(psi) << " ";
    for(int j = 0; j<N; j++){ //save local energy values
      enerfile << LocalEnergy[j] << " ";
    }
    enerfile << std::endl;

    printfln("Iteration %d,energy = %0.3g, max link dim is %d",n+1,energy,maxLinkDim(psi));

  }

  //calculate properties of ground state of final H
  auto [finalEnergy, psiF] = dmrg(toMPO(ampo),MPS(state),sweeps,{"Silent=",true});
  LocalEnergy = localEnergy(Nx, Ny, sites, psiF, LED, LEDyPBC);
  enerfile << 0.0 << " " << finalEnergy << " " << maxLinkDim(psiF) << " ";
  for(int j = 0; j<N; j++){ //save local energy values
    enerfile << LocalEnergy[j] << " ";
  }
    enerfile << std::endl;

  enerfile.close();

  return 0;
  }

// calculates local energy for 2D MPS using gates
std::vector<double> localEnergy(int Nx, int Ny, SiteSet sites, MPS psi, std::vector<std::vector<ITensor>> LED, std::vector<ITensor> LEDyPBC){
  std::vector<double> LocalEnergy(length(psi),0.0);
  for(int i=1; i<=Nx; i++){
    for(int j=1; j<=Ny; j++){ 
      int index = (i-1)*Ny + j; //this order to make orthogonality center easier to compute
      ITensor ket;
      if(j==Ny){ //y-periodic boundary equations with swap gates
        psi.position(index-Ny+1);
        for(int n=1; n<Ny-1; n++){
          int b = index-Ny+n; //define gate bond
          auto g = BondGate(sites,b,b+1);
          auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
          psi.position(g.i2()); //orthogonality center moves to the right
        }
        
        ket = psi(index-1)*psi(index);
        LocalEnergy[index-1] = eltC( dag(prime(ket,"Site")) * LEDyPBC[i-1] * ket).real();

        //restore the state to the original MPS
        for(int n=Ny-1; n>1; n--){
          int b = index-Ny+n;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromright); //svd from the right
          psi.position(g.i1()); //move orthogonality center to the left  
        }
      }// y-periodic

      //original nearest-neighbour code
      else{ 
        psi.position(index);
        ket = psi(index)*psi(index+1);
        LocalEnergy[index-1] = eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();
      } //nearest-neighbour

      if(i<Nx){
        psi.position(index+Ny); //site to bring to index+1  

        // bring index+Ny to position index+1
        for(int n=Ny; n>1; n--){
          int b = index+n;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
          psi.position(g.i1()); //orthogonality center moves to the left
        }

        ket = psi(index)*psi(index+1);
        LocalEnergy[index-1] += eltC( dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();

        // bring index+1 back to position index+Ny
        for(int n=1; n<Ny; n++){
          int b = index+n;
          auto g = BondGate(sites,b,b+1);
          auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
          psi.position(g.i2()); //orthogonality center moves to the right
        }
      }//long-range interaction
    }// for j
  }// for i
  psi.orthogonalize({"Cutoff=",1E-12}); //restore psi bond dimensions
  return LocalEnergy;
}//localEnergy
