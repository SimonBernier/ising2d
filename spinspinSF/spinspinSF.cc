#include "itensor/all.h"

using namespace itensor;

//time stepping by dt with MPO method
void timeStepMPO(AutoMPO, MPS&, double, Args);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, int, double, SiteSet, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);

int main(int argc, char *argv[])
  {
  int Ly = 4;
  int Lx = 16;
  float h = 4.0;
  int Nt = 20;
  float dt = 0.1;
  int method = 0; //0 for MPO, 1 for TEBD

  if(argc > 6)
    method = std::stoi(argv[6]);
  if(argc > 5)
    dt = std::stof(argv[5]);
  if(argc > 4)
    Nt = std::stoi(argv[4]);
  if(argc > 3)
    h = std::stof(argv[3]);
  if(argc > 2)
    Lx = std::stoi(argv[2]);
  if(argc > 1)
    Ly = std::stoi(argv[1]);

  // write results to file
  char schar[64];
  if(method == 0){
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_dt_%0.3g_Nt_%d_mpo.dat",Ly,Lx,h,dt,Nt); 
  } else if(method == 1){
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_dt_%0.3g_Nt_%d_tebd2.dat",Ly,Lx,h,dt,Nt);
  } else if(method == 2){
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_dt_%0.3g_Nt_%d_tebd4.dat",Ly,Lx,h,dt,Nt);
  } else{
    printfln("Not a valid method");
    return 0;
  }
  std::string s1(schar);
  std::ofstream dataFile;
  dataFile.open(s1); // opens the file
  if( !dataFile ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
  }
  //make header for t=0 calculations
  dataFile << "t=0" << " " << "enPsi" << " " << "MaxDimPsi" << " " << "enPhi" << " " << "MaxDimPhi" << " " 
           << "Sz(x,y)" << " " << "Sz(x,y)Sz(Lx/2,1,0)" << " " << std::endl;

  auto L = Ly * Lx;
  auto sites = SpinHalf(L,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

  // autompo hamiltonian
  for(auto j : lattice){
    ampo += -4, "Sz", j.s1, "Sz", j.s2;
  }
  for(auto j : range1(L)){
    ampo += -2.0*h, "Sx", j;
  }
  auto H = toMPO(ampo);

  //initial state
  auto state = InitState(sites); 
  for(auto j : range1(L)){
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
  }

  // 2d ising model parameters
  auto sweeps = Sweeps(10);
  sweeps.maxdim() = 20, 50, 100, 100, 200, 200, 512;
  sweeps.cutoff() = 1E-10;
  
  //make 2D vector of ITensor for local energy operators
  //long-range interactions have the same structure as nearest-neighbour when we use swap gates
  std::vector<std::vector<ITensor>> LED(Lx, std::vector<ITensor>(Ly-1));
  std::vector<std::vector<ITensor>> LED_LR(Lx-1, std::vector<ITensor>(Ly));
  std::vector<ITensor> LEDyPBC(Lx);
  // make local energy tensors
  for(int i=1; i<=Lx; i++){
    for(int j=1; j<=Ly; j++){
      int index = (i-1)*Ly+j;
      //MPS long-range
      if(i<Lx && j==1){
        for(int m = 0; m<Ly; m++){
          LED_LR[i-1][m] = -4.0*sites.op("Sz",index+2*m)*sites.op("Sz",index+2*m+1);
        }
      }
      //y-periodic boundary equations
      if(j==Ly){
        // site index-Ly+1 is moved to site index-1 with swap gates
        LEDyPBC[i-1] = -4.0*sites.op("Sz",index-1)*sites.op("Sz",index);
        LEDyPBC[i-1] += -2.0*h*sites.op("Id",index-1)*sites.op("Sx",index);
      }
      // MPS nearest-neighbour
      if(j<Ly){
        LED[i-1][j-1] = -4.0*sites.op("Sz",index)*sites.op("Sz",index+1);
        LED[i-1][j-1] += -2.0*h*sites.op("Sx",index)*sites.op("Id",index+1);
      }
    }
  }

  //make vector of Sz MPOs
  std::vector<MPO> Sz(L);
  for(auto j : range1(L)){
    auto ampoSz = AutoMPO(sites);
    for(auto i : range1(L)){
      if(i==j){
        ampoSz += 2.0, "Sz", i; //Sz only at site j
      }
      else{
        ampoSz += "Id", i; //identities everywhere else
      }
    }
    Sz[j-1] = toMPO(ampoSz);
  }

  // calculate ground state
  auto [en_psi, psi] = dmrg(H,MPS(state),sweeps,{"Silent=",true});

  //calculate <Sz(x,y)>
  std::vector<double> aveSz(L,0.0);
  for(auto j : range1(L)){
    psi.position(j);
    aveSz[j-1] = elt( dag(prime(psi(j),"Site")) * 2.0*sites.op("Sz",j) * psi(j) );
  }

  // make |phi> = Sz|psi>
  int loc = (Lx/2)*Ly+1; //centered in x and on lower row in y
  psi.position(loc);
  auto newA = 2.0*sites.op("Sz",loc) * psi(loc);
  newA.noPrime();
  auto phi = psi;
  phi.set(loc, newA);
  auto en_phi = inner(phi,H,phi);

  // store spin-spin correlation function
  std::vector<Complex> szsz(L,0.0); // local energy density vector
  //calculate the spin-spin correlation using MPS * MPO * MPS methods
  for(auto j : range1(L)){
    szsz[j-1] = innerC(psi, Sz[j-1], phi);
  }

  printfln("\nIteration %d, time = %0.2f; phi energy = %0.f, max link dim is %d",0,0, en_phi,maxLinkDim(phi));
  // store to file
  dataFile << 0.0 << " " << en_psi << " " << maxLinkDim(psi) << " " << en_phi << " " << maxLinkDim(phi) << " ";
  for(int j = 0; j<L; j++){ //save local energy values
    dataFile << aveSz[j] << " ";
  }
  for(int j = 0; j<L; j++){ //save local energy values
    dataFile << real(szsz[j]) << " " << imag(szsz[j]) << " ";
  }
  dataFile << std::endl;

  dataFile << "tval" << " " << "enPhi" << " " << "MaxDimPhi" << " " << "Sz(x,y,t)Sz(Lx/2,1,0)" << " " << std::endl;

  // time evolution parameters
  double tval = 0.0; //time

  //args for time evolution methods
  Args args;
  std::vector<BondGate> gates; //only make the gates vector if using TEBD
  if(method==0){
    printfln("Starting MPO based time evolution");
    args = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",512);
  } else if(method==1){
    printfln("Starting second order TEBD");
    args = Args("Cutoff=",1E-10,"MaxDim=",512);
    //Create a std::vector (dynamically sizeable array) to hold the Trotter gates
    gates = makeGates(Lx, Ly, dt, sites, LED, LEDyPBC, LED_LR);
  } else if(method==2){
    printfln("Starting fourth order TEBD");
    args = Args("Cutoff=",1E-10,"MaxDim=",512);
    //Create a std::vector (dynamically sizeable array) to hold the Trotter gates
    Real delta1 =  0.414490771794376*dt;
    Real delta2 = -0.657963087177503*dt;
    auto gatesdelta1 = makeGates(Lx, Ly, delta1, sites, LED, LEDyPBC, LED_LR);
    auto gatesdelta2 = makeGates(Lx, Ly, delta2, sites, LED, LEDyPBC, LED_LR);
    gates = gatesdelta1;
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
    gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////// time evolve //////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  for(int n=1; n<=Nt; n++){
    tval += dt; //update time vector
    // do a time step
    if(method == 0){
      // MPO time step, overwriting phi when done
      timeStepMPO(ampo, phi, dt, args);
    }
    else if(method == 1){
      //Time evolve, orthogonalizing and overwriting phi when done
      gateTEvol(gates,dt,dt,phi,{args,"Verbose=",false});
      phi.orthogonalize(args);
    }
    else if(method == 2){
      //4th order TEBD
      gateTEvol(gates,dt,dt,phi,{args,"Verbose=",false});
      phi.orthogonalize(args);
    }

    en_phi = innerC(phi, H, phi).real();

    //calculate the spin-spin correlation using MPS * MPO * MPS methods
    for(auto j : range1(L)){
      szsz[j-1] = innerC(exp(1_i*en_psi*tval)*psi, Sz[j-1], phi);
    }

    //write to file
    dataFile << tval << " " << en_phi << " " << maxLinkDim(phi) << " ";
    for(int j = 0; j<L; j++){ //save local energy values
      dataFile << real(szsz[j]) << " " << imag(szsz[j]) << " ";
    }
    dataFile << std::endl;

    printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(phi));

  }

  dataFile.close();

  return 0;
  }

// do one MPO based time evolution step
void timeStepMPO(AutoMPO ampo, MPS& psi, double dt, Args args)
  {
    //time evolution operators
    auto expH1 = toExpH(ampo, 0.5*dt*(Cplx_i+1)); //time evolve by -0.5*(i+1)*dt
    auto expH2 = toExpH(ampo, 0.5*dt*(Cplx_i-1)); //time evolve by -0.5*(i-1)*dt

    //Fit method for MPO*MPS
    psi = applyMPO(expH1,psi,args);
    psi.noPrime(); //need to do this after each to take care of prime levels
    psi = applyMPO(expH2,psi,args);
    psi.noPrime().normalize(); //need to do this after each to take care of prime levels

  }

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int Lx, int Ly, double dt, SiteSet sites, 
                                std::vector<std::vector<ITensor>> LED,
                                std::vector<ITensor> LEDyPBC,
                                std::vector<std::vector<ITensor>> LED_LR)
  {
  std::vector<BondGate> gates; 
  //Create the gates exp(-i*tstep/2*hterm)
  for(int i=1; i<=Lx; i++){
    for(int j=1; j<=Ly; j++){ 
      int index = (i-1)*Ly + j; //MPS site index

      if(j==1){ //y-periodic boundary equations with swap gates
        for(int n=0; n<Ly-2; n++){ //swap from index-Ly+1 to index-1
          int b = index+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
        }
        auto hterm = LEDyPBC[i-1];
        auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        //restore the state to the original MPS
        for(int n=Ly-2; n>0; n--){
          int b = index+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }
      }// y-periodic

        //original nearest-neighbour code
      if(j<Ly){
        auto hterm = LED[i-1][j-1];
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
      } //nearest-neighbour

      // long-range interaction
      if(i<Lx && j==1){ // bring index+Ly to position index+1
        for(int m=0; m<=Ly-2; m++){
          for(int n=Ly; n>1+m; n--){
            int b = index + n + m;
            auto swapGate = BondGate(sites,b-1,b);
            gates.push_back(swapGate);
          }
        }

        for(int m = 0; m<Ly; m++){
          auto hterm = LED_LR[i-1][m];
          auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
          gates.push_back(g);
        }

        // bring index+1 back to position index+Ly
        for(int m=Ly-2; m>=0; m--){
          for(int n=1+m; n<Ly; n++){
            int b = index + n + m;
            auto swapGate = BondGate(sites,b,b+1);
            gates.push_back(swapGate);
          }
        }
      }//long-range interaction
    }// for j
  }// for i

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int i=Lx; i>=1; i--){
    for(int j=Ly; j>=1; j--){ 
      int index = (i-1)*Ly + j; //MPS site index

      // long-range interaction
      if(i<Lx && j==1){ // bring index+Ly to position index+1
        for(int m=0; m<=Ly-2; m++){
          for(int n=Ly; n>1+m; n--){
            int b = index + n + m;
            auto swapGate = BondGate(sites,b-1,b);
            gates.push_back(swapGate);
          }
        }
        
        for(int m = Ly-1; m>=0; m--){
          auto hterm = LED_LR[i-1][m];
          auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
          gates.push_back(g);
        }

        // bring index+1 back to position index+Ly
        for(int m=Ly-2; m>=0; m--){
          for(int n=1+m; n<Ly; n++){
            int b = index + n + m;
            auto swapGate = BondGate(sites,b,b+1);
            gates.push_back(swapGate);
          }
        }
      }//long-range interaction

      //original nearest-neighbour code
      if(j<Ly){
        auto hterm = LED[i-1][j-1];
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
      } //nearest-neighbour

      //y-periodic boundary equations with swap gates
      if(j==1){
        for(int n=0; n<Ly-2; n++){ //swap from index-Ly+1 to index-1
          int b = index+n;
          auto swapGate = BondGate(sites,b,b+1);
          gates.push_back(swapGate);
        }
        auto hterm = LEDyPBC[i-1];
        auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);

        //restore the state to the original MPS
        for(int n=Ly-2; n>0; n--){
          int b = index+n;
          auto swapGate = BondGate(sites,b-1,b);
          gates.push_back(swapGate);
        }
      }// y-periodic
    }// for j
  }// for i

  return gates;
  
}// makeGates