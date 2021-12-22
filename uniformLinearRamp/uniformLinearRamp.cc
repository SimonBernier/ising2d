#include "itensor/all.h"

using namespace itensor;

//function definition for calculation of local energy
std::vector<double> localEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);
//time stepping by dt with MPO method
void timeStepMPO(AutoMPO, MPS&, double, Args);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, int, double, double, SiteSet, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);
// function to update time dependent parts of the gates vector
void updateGates(int, int, double, double, SiteSet, std::vector<BondGate>&, std::vector<std::vector<ITensor>>,
                    std::vector<ITensor>, std::vector<std::vector<ITensor>>);


int main(int argc, char *argv[])
  {
  int Ly = 4;
  int Lx = 16;
  float h = 4.0;
  int Nt = 20;
  float dt = 0.1;
  float tau = 1.0;
  int method = 0; //0 for MPO, 1 for TEBD

  if(argc > 7)
    method = std::stoi(argv[7]);
  if(argc > 6)
    tau = std::stof(argv[6]);
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
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_tau_%0.3g_dt_%0.3g_mpo.dat",Ly,Lx,h,tau,dt); 
  } else if(method == 1){
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_tau_%0.3g_dt_%0.3g_tebd.dat",Ly,Lx,h,tau,dt);
  } else{
    printfln("Not a valid method");
    return 0;
  }
  std::string s1(schar);
  std::ofstream enerfile;
  enerfile.open(s1); // opens the file
  if( !enerfile ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
  }
  //make header
  enerfile << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;

  auto L = Ly * Lx;
  auto sites = SpinHalf(L,{"ConserveQNs=",false});

  auto ampo = AutoMPO(sites);
  auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

  // autompo hamiltonian
  for(auto j : lattice){
      ampo += -4, "Sz", j.s1, "Sz", j.s2;
  }
  //final Hamiltonian
  double hF = 2.6; //final field (gapless ground state)
  for(auto j : range1(L)){
    ampo += -2.0*hF, "Sx", j;
  }
  auto Hf = toMPO(ampo);

  //initial state
  auto state = InitState(sites); 
  for(auto j : range1(L)){
      state.set(j, (j % 2 == 1 ? "Up" : "Dn"));
  }

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200, 400;
  sweeps.cutoff() = 1E-10;

  // calculate initial local energy density
  std::vector<double> LocalEnergy(L,0.0); // local energy density vector
  
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
        LEDyPBC[i-1] += -2.0*hF*sites.op("Id",index-1)*sites.op("Sx",index);
      }
      // MPS nearest-neighbour
      if(j<Ly){
        LED[i-1][j-1] = -4.0*sites.op("Sz",index)*sites.op("Sz",index+1);
        LED[i-1][j-1] += -2.0*hF*sites.op("Sx",index)*sites.op("Id",index+1);
      }
    }
  }

  //// EVENTUALLY, THIS PART OF THE CODE MUST BE COMMENTED OUT
  //// CALCULATING THE CRITICAL GROUND STATE IS DIFFICULT FOR
  //// LARGE SYSTEM SIZES. KEEP IT FOR NOW FOR TESTING
  // calculate ground state of final H
  auto [finalEnergy, psiF] = dmrg(Hf,MPS(state),sweeps,{"Silent=",true});
  // calculate local energy
  LocalEnergy = localEnergy(Lx, Ly, sites, psiF, LED, LEDyPBC, LED_LR);
  // store to file
  enerfile << 0.0 << " " << finalEnergy << " " << maxLinkDim(psiF) << " ";
  for(int j = 0; j<L; j++){ //save local energy values
  enerfile << LocalEnergy[j] << " ";
  }
  enerfile << std::endl;

  // revert to H at t=0 (gapped ground state)
  for(auto j : range1(L)){
    ampo += -2.0*(h-hF), "Sx", j;
  }

  //DMRG to find ground state at t=0
  auto [energy,psi] = dmrg(toMPO(ampo),MPS(state),sweeps,{"Silent=",true});
  energy = inner(psi,Hf,psi); //energy <psi(0)|H_final|psi(0)>
  // calculate local energy density
  LocalEnergy = localEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);
  //store to file
  enerfile << 0.0 << " " << energy << " " << maxLinkDim(psi) << " "; //print to file
  for(int j = 0; j<L; j++){ //save local energy values
    enerfile << LocalEnergy[j] << " ";
  }
  enerfile << std::endl;

  // time evolution parameters
  double tval = 0.0; //time
  auto hval = std::vector<double>(Nt+1,0.0);
  auto diff = std::vector<double>(Nt+1,0.0);
  int Nstep = tau/dt;
  double step = (h-hF)/Nstep;
  hval[0] = h;
  diff[0] = 0.0;
  for(int j=1; j<=Nt; j++){ //make linear ramp
    if(j<=Nstep){
      hval[j] = h-step*j; //printf("h = %0.2f, ", hval[j]);
    } else{
      hval[j] = hval[j-1]; //printf("h = %0.2f, ", hval[j]);
    }
    diff[j] = hval[j]-hval[j-1]; //printfln("diff = %0.2f, ", diff[j]);
  }
  //args for time evolution methods
  Args args;
  std::vector<BondGate> gates; //only make the gates vector if using TEBD
  if(method==0){
    printfln("Starting MPO based time evolution");
    args = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",3000);
  } else if(method==1){
    printfln("Starting TEBD");
    args = Args("Cutoff=",1E-10,"MaxDim=",3000);
    //Create a std::vector (dynamically sizeable array) to hold the Trotter gates
    gates = makeGates(Lx, Ly, hval[0], dt, sites, LED, LEDyPBC, LED_LR);
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////// time evolve //////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  for(int n=1; n<=Nt; n++){
    // do a time step
    if(method == 0){
      // update autoMPO
      for(auto j : range1(L)){
        ampo += -2.0*diff[n], "Sx", j;
      }
      // MPO time step, overwriting psi when done
      timeStepMPO(ampo, psi, dt, args);
    }
    else if(method == 1){
      //Create a std::vector (dynamically sizeable array) to hold the Trotter gates
      updateGates(Lx, Ly, hval[n]-hF, dt, sites, gates, LED, LEDyPBC, LED_LR);
      //Time evolve, orthogonalizing and overwriting psi when done
      gateTEvol(gates,dt,dt,psi,{args,"Verbose=",false,"Normalize=",true});
      psi.orthogonalize(args);
    }

    // calculate energy <psi(t)|H_final|psi(t)>
    energy = innerC(psi,Hf,psi).real();
    // calculate local energy density <psi(t)|H_final(x,y)|psi(t)>
    LocalEnergy = localEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);

    //write to file
    tval += dt; //update time vector
    enerfile << tval << " " << energy << " " << maxLinkDim(psi) << " ";
    for(int j = 0; j<L; j++){ //save local energy values
      enerfile << LocalEnergy[j] << " ";
    }
    enerfile << std::endl;

    printfln("\nIteration %d, time = %0.2f; energy = %0.3g, max link dim is %d",n,tval,energy,maxLinkDim(psi));

  }

  enerfile.close();

  return 0;
  }

// calculates local energy for 2D MPS using gates
std::vector<double> localEnergy(int Lx, int Ly, SiteSet sites, MPS psi,
                                std::vector<std::vector<ITensor>> LED, 
                                std::vector<ITensor> LEDyPBC,
                                std::vector<std::vector<ITensor>> LED_LR)
  {  
  std::vector<double> LocalEnergy(length(psi),0.0);
  for(int i=1; i<=Lx; i++){
    for(int j=1; j<=Ly; j++){ 
      int index = (i-1)*Ly + j; //this order to make orthogonality center easier to compute
      ITensor ket;
      if(j==1){ //y-periodic boundary equations with swap gates
        psi.position(index);
        for(int n=0; n<Ly-2; n++){
          int b = index+n;//define gate bond
          auto g = BondGate(sites,b,b+1);
          auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
          psi.position(g.i2()); //orthogonality center moves to the right
        }
        
        ket = psi(index+Ly-2)*psi(index+Ly-1);
        LocalEnergy[index+Ly-2] += eltC( dag(prime(ket,"Site")) * LEDyPBC[i-1] * ket).real();

        //restore the state to the original MPS
        for(int n=Ly-2; n>0; n--){
          int b = index+n;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromright); //svd from the right
          psi.position(g.i1()); //move orthogonality center to the left  
        }
      }// y-periodic

      //original nearest-neighbour code
      if(j<Ly){ 
        psi.position(index);
        ket = psi(index)*psi(index+1);
        LocalEnergy[index-1] += eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();
      } //nearest-neighbour

      //smart ordering of gates for sites i*Ly+1 to i*Ly+(Ly-1) with sites (i+1)*Ly+1 to (i+1)*Ly+(Ly-1)
      if(i<Lx && j==1){  
        // bring index+Ly to position index+1
        for(int m=0; m<=Ly-2; m++){
          psi.position(index+Ly+m);
          for(int n=Ly; n>1+m; n--){
          int b = index + n + m;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
          psi.position(g.i1()); //orthogonality center moves to the left
          }
        }
        
        for(int m = 0; m<Ly; m++){
          psi.position(index+2*m);
          ket = psi(index+2*m)*psi(index+2*m+1);
          LocalEnergy[index-1+m] += eltC( dag(prime(ket,"Site")) * LED_LR[i-1][m] * ket).real();
        }

        // bring index+1 back to position index+Ly
        for(int m=Ly-2; m>=0; m--){
          psi.position(index+1+2*m);
          for(int n=1+m; n<Ly; n++){
            int b = index + n + m;
            auto g = BondGate(sites,b,b+1);
            auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi.position(g.i2()); //orthogonality center moves to the right
          }
        }
      }//long-range interaction
    }// for j
  }// for i

  return LocalEnergy;

}//localEnergy

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
std::vector<BondGate> makeGates(int Lx, int Ly, double hDiff, double dt, SiteSet sites, 
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
        hterm += -2.0*hDiff*op(sites,"Id",index+Ly-2)*op(sites,"Sx",index+Ly-1);
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
        hterm += -2.0*hDiff*op(sites,"Sx",index)*op(sites,"Id",index+1);
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
        hterm += -2.0*hDiff*op(sites,"Sx",index)*op(sites,"Id",index+1);
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
        hterm += -2.0*hDiff*op(sites,"Id",index+Ly-2)*op(sites,"Sx",index+Ly-1);
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

//update time dependent parts of the gates
void updateGates(int Lx, int Ly, double hDiff, double dt, SiteSet sites, std::vector<BondGate>& gates, 
                  std::vector<std::vector<ITensor>> LED, std::vector<ITensor> LEDyPBC, std::vector<std::vector<ITensor>> LED_LR)
  {
  int indG = 0; //location of time dependent parts of gates
  //Create the gates exp(-i*tstep/2*hterm)
  for(int i=1; i<=Lx; i++){
    for(int j=1; j<=Ly; j++){ 
      int index = (i-1)*Ly + j; //MPS site index

      if(j==1){ //y-periodic boundary equations with swap gates
        indG += Ly-2; //swap gates
        auto hterm = LEDyPBC[i-1];
        hterm += -2.0*hDiff*op(sites,"Id",index+Ly-2)*op(sites,"Sx",index+Ly-1);
        auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
        gates[indG] = g;
        indG++; //t-dep gate
        indG += Ly-2; //swapgates
      }// y-periodic

      //original nearest-neighbour code
      if(j<Ly){
        auto hterm = LED[i-1][j-1];
        hterm += -2.0*hDiff*op(sites,"Sx",index)*op(sites,"Id",index+1);
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates[indG] = g;
        indG++;
      } //nearest-neighbour

      // long-range interaction
      if(i<Lx && j==1){ // bring index+Ly to position index+1
        indG += Ly*(Ly-1)/2; //swap gates
        for(int m = 0; m<Ly; m++){
          auto hterm = LED_LR[i-1][m];
          auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
          gates[indG] = g;
          indG++;
        }
        indG += Ly*(Ly-1)/2; //swap gates
      }//long-range interaction
    }// for j
  }// for i

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int i=Lx; i>=1; i--){
    for(int j=Ly; j>=1; j--){ 
      int index = (i-1)*Ly + j; //MPS site index

      // long-range interaction
      if(i<Lx && j==1){ // bring index+Ly to position index+1
        indG += Ly*(Ly-1)/2; //swap gates       
        for(int m = Ly-1; m>=0; m--){
          auto hterm = LED_LR[i-1][m];
          auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
          gates[indG] = g;
          indG++;
        }
        indG += Ly*(Ly-1)/2; //swap gates
      }//long-range interaction

      //original nearest-neighbour code
      if(j<Ly){
        auto hterm = LED[i-1][j-1];
        hterm += -2.0*hDiff*op(sites,"Sx",index)*op(sites,"Id",index+1);
        auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
        gates[indG] = g;
        indG++;
      } //nearest-neighbour

      //y-periodic boundary equations with swap gates
      if(j==1){
        indG += Ly-2;
        auto hterm = LEDyPBC[i-1];
        hterm += -2.0*hDiff*op(sites,"Id",index+Ly-2)*op(sites,"Sx",index+Ly-1);
        auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
        gates[indG] = g;
        indG++;
        indG += Ly-2;
      }// y-periodic
    }// for j
  }// for i

  return;

  }