#include "itensor/all.h"

using namespace itensor;

//function definition for calculation of local energy
std::vector<double> localEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>, std::vector<ITensor>);
//function which makes the bond gates to pass to gateTEvol
std::vector<BondGate> makeGates(int, int, double, double, SiteSet, std::vector<std::vector<ITensor>>, std::vector<ITensor>);


int main(int argc, char *argv[]){

  int Nx = 16;
  int Ny = 4;
  double h = 4.0;
  int Nt = 20;
  double dt = 0.1;
  double tau = 1.0;
  if(argc > 6)
    tau = std::stod(argv[6]);
  if(argc > 5)
    dt = std::stod(argv[5]);
  if(argc > 4)
    Nt = std::stoi(argv[4]);
  if(argc > 3)
    h = std::stod(argv[3]);
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
  auto args = Args("Cutoff=",1E-10,"Verbose=",false);

  //time evolution fields h and diff
  double tval = 0.0; //time
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
  for(int n=1; n<=Nt; n++){
    //
    // update autoMPO for energy calculation
    for(auto j : range1(N)){
      ampo += -2.0*diff[n], "Sx", j;
    }

    //Create a std::vector (dynamically sizeable array)
    //to hold the Trotter gates
    std::vector<BondGate> gates = makeGates(Nx, Ny, hval[n], dt, sites, LED, LEDyPBC);

    //Time evolve, overwriting psi when done
    gateTEvol(gates,dt,dt,psi,args); //evolve only once since gates are time dependent
    psi.orthogonalize({"Cutoff=",1E-12}); //restore psi bond dimensions
    tval += dt; //update time

    //Calculate ground state energy
    auto energy = innerC(psi,toMPO(ampo),psi).real();

    // calculate local energy density
    LocalEnergy = localEnergy(Nx, Ny, sites, psi, LED, LEDyPBC);

    //write to file
    enerfile << tval << " " << energy << " " << maxLinkDim(psi) << " ";
    for(int j = 0; j<N; j++){ //save local energy values
      enerfile << LocalEnergy[j] << " ";
    }
    enerfile << std::endl;

    printfln("\nIteration %d, time = %0.2f; energy = %0.3g, max link dim is %d",n,tval,energy,maxLinkDim(psi));

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

}//main

// calculates local energy for 2D MPS using gates
std::vector<double> localEnergy(int Nx, int Ny, SiteSet sites, MPS psi, std::vector<std::vector<ITensor>> LED, std::vector<ITensor> LEDyPBC)
  {  
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

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int Nx, int Ny, double h, double dt, SiteSet sites, std::vector<std::vector<ITensor>> LED, std::vector<ITensor> LEDyPBC)
  {
  std::vector<BondGate> gates; 
  //Create the gates exp(-i*tstep/2*hterm)
  for(int i=1; i<=Nx; i++){
    for(int j=1; j<=Ny; j++){ 
      int index = (i-1)*Ny + j; //MPS site index
      if(j==Ny){ //y-periodic boundary equations with swap gates
        for(int n=1; n<Ny-1; n++){ //swap from index-Ny+1 to index-1
          int b = index-Ny+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
         }
        auto hterm = LEDyPBC[i-1];
        hterm += -2.0*h*op(sites,"Id",index-1)*op(sites,"Sx",index);
        auto g = BondGate(sites,index-1,index,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        //restore the state to the original MPS
        for(int n=Ny-1; n>1; n--){
          int b = index-Ny+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }
      }// y-periodic

        //original nearest-neighbour code
        else{
        auto hterm = LED[i-1][j-1];
        hterm += -2.0*h*op(sites,"Sx",index)*op(sites,"Id",index+1);
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
      } //nearest-neighbour

        // long-range interaction
      if(i<Nx){ // bring index+Ny to position index+1
        for(int n=Ny; n>1; n--){
          int b = index+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }

        auto hterm = LED[i-1][j-1];
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        // bring index+1 back to position index+Ny
        for(int n=1; n<Ny; n++){
          int b = index+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
        }
      }//long-range interaction
    }// for j
  }// for i

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int i=Nx; i>=1; i--){
    for(int j=Ny; j>=1; j--){ 
      int index = (i-1)*Ny + j; //MPS site index

      // long-range interaction
      if(i<Nx){ // bring index+Ny to position index+1
        for(int n=Ny; n>1; n--){
          int b = index+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }

        auto hterm = LED[i-1][j-1];
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        // bring index+1 back to position index+Ny
        for(int n=1; n<Ny; n++){
          int b = index+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
        }
      }//long-range interaction

      //y-periodic boundary equations with swap gates
      if(j==Ny){
        for(int n=1; n<Ny-1; n++){ //swap from index-Ny+1 to index-1
          int b = index-Ny+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
        }
        auto hterm = LEDyPBC[i-1];
        hterm += -2.0*h*op(sites,"Id",index-1)*op(sites,"Sx",index);
        auto g = BondGate(sites,index-1,index,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        //restore the state to the original MPS
        for(int n=Ny-1; n>1; n--){
          int b = index-Ny+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }
      }// y-periodic
      //original nearest-neighbour code
      else{
        auto hterm = LED[i-1][j-1];
        hterm += -2.0*h*op(sites,"Sx",index)*op(sites,"Id",index+1);
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
      } //nearest-neighbour
    }// for j
  }// for i

  return gates;
  
}// makeGates