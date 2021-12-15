#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
  {
  int Nx = 16;
  int Ny = 4;
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
  auto [energy,psi0] = dmrg(H,MPS(state),sweeps,{"Silent=",true});
  
  // calculate local energy density
  for(int i=1; i<=Nx; i++){
    for(int j=1; j<=Ny; j++){ 
      int index = (i-1)*Ny + j; //this order to make orthogonality center easier to compute
      ITensor ket;
      if(j==Ny){ //y-periodic boundary equations with swap gates
        psi0.position(index-Ny+1);
        for(int n=1; n<Ny-1; n++){
          int b = index-Ny+n; //define gate bond
          auto g = BondGate(sites,b,b+1);
          auto AA = psi0(b)*psi0(b+1)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi0.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
          psi0.position(g.i2()); //orthogonality center moves to the right
        }
        
        ket = psi0(index-1)*psi0(index);
        auto energy = elt( dag(prime(ket,"Site")) * LEDyPBC[i-1] * ket);
        LocalEnergy[index-1] = energy;

        //restore the state to the original MPS
        for(int n=Ny-1; n>1; n--){
          int b = index-Ny+n;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi0(b-1)*psi0(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0");
          psi0.svdBond(g.i1(), AA, Fromright); //svd from the right
          psi0.position(g.i1()); //move orthogonality center to the left  
        }
      }// y-periodic

      //original nearest-neighbour code
      else{ 
        psi0.position(index);
        ket = psi0(index)*psi0(index+1);
        auto energy = elt(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket);
        LocalEnergy[index-1] = energy;
      } //nearest-neighbour

      if(i<Nx){
        psi0.position(index+Ny); //site to bring to index+1  

        // bring index+Ny to position index+1
        for(int n=Ny; n>1; n--){
          int b = index+n;
          auto g = BondGate(sites,b-1,b);
          auto AA = psi0(b-1)*psi0(b)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi0.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
          psi0.position(g.i1()); //orthogonality center moves to the left
        }

        ket = psi0(index)*psi0(index+1);
        auto energy = elt(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket);
        LocalEnergy[index-1] += energy;

        // bring index+1 back to position index+Ny
        for(int n=1; n<Ny; n++){
          int b = index+n;
          auto g = BondGate(sites,b,b+1);
          auto AA = psi0(b)*psi0(b+1)*g.gate(); //contract over bond b
          AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
          psi0.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
          psi0.position(g.i2()); //orthogonality center moves to the right
        }
      }//long-range interaction
    }// for j
  }// for i
  psi0.orthogonalize({"Cutoff=",1E-12}); //restore psi0 bond dimensions

  // create vectors of time
  int Nt = 20;
  double dt = 0.1;
  std::vector<double> tval(Nt+1, 0.0);
  
  enerfile1 << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;
  enerfile2 << "tval" << " " << "energy" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;
  enerfile1 << tval[0] << " " << energy << " " << maxLinkDim(psi0) << " "; //print to file
  enerfile2 << tval[0] << " " << energy << " " << maxLinkDim(psi0) << " "; //print to file
  for(int j = 0; j<N; j++){ //save local energy values
    enerfile1 << LocalEnergy[j] << " ";
    enerfile2 << LocalEnergy[j] << " ";
  }
  enerfile1 << std::endl;
  enerfile2 << std::endl;


  // time evolution parameters 
  auto args_DM = Args("Method=","DensityMatrix","Cutoff=",1E-10,"MaxDim=",3000);
  auto args_Fit = Args("Method=","Fit","Cutoff=",1E-10,"MaxDim=",3000);

  //time evolution fields h and diff
  auto hval = std::vector<double>(Nt,0.0);
  auto diff = std::vector<double>(Nt,0.0);
  double step = 2.0*(h-2.0)/Nt;
  hval[0] = h-step;
  diff[0] = -step;
  for(int j=1; j<Nt; j++){ //make linear ramp
    if(j<=Nt/2){
      hval[j] = h-step*j;
    } else{
      hval[j] = hval[j-1];
    }
    diff[j] = hval[j]-hval[j-1];
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
    tval[i+1]=tval[i]+dt; //update time vector

    // check DensityMatrix method for MPO*MPS
    psi_DM = applyMPO(expH1,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_DM = applyMPO(expH2,psi_DM,args_DM);
    psi_DM.noPrime().normalize(); //need to do this after each to take care of prime levels
    auto energy_DM = innerC(psi_DM,H,psi_DM).real();

    // calculate local energy density
    for(int i=1; i<=Nx; i++){
      for(int j=1; j<=Ny; j++){ 
        int index = (i-1)*Ny + j; //this order to make orthogonality center easier to compute
        ITensor ket;
        if(j==Ny){ //y-periodic boundary equations with swap gates
          psi_DM.position(index-Ny+1);
          for(int n=1; n<Ny-1; n++){
            int b = index-Ny+n; //define gate bond
            auto g = BondGate(sites,b,b+1);
            auto AA = psi_DM(b)*psi_DM(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_DM.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi_DM.position(g.i2()); //orthogonality center moves to the right
          }
          
          ket = psi_DM(index-1)*psi_DM(index);
          LocalEnergy[index-1] = eltC( dag(prime(ket,"Site")) * LEDyPBC[i-1] * ket).real();;

          //restore the state to the original MPS
          for(int n=Ny-1; n>1; n--){
            int b = index-Ny+n;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi_DM(b-1)*psi_DM(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0");
            psi_DM.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi_DM.position(g.i1()); //move orthogonality center to the left  
          }
        }// y-periodic

        //original nearest-neighbour code
        else{ 
          psi_DM.position(index);
          ket = psi_DM(index)*psi_DM(index+1);
          LocalEnergy[index-1] = eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();
        } //nearest-neighbour

        if(i<Nx){
          psi_DM.position(index+Ny); //site to bring to index+1  

          // bring index+Ny to position index+1
          for(int n=Ny; n>1; n--){
            int b = index+n;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi_DM(b-1)*psi_DM(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_DM.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi_DM.position(g.i1()); //orthogonality center moves to the left
          }

          ket = psi_DM(index)*psi_DM(index+1);
          LocalEnergy[index-1] += eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();

          // bring index+1 back to position index+Ny
          for(int n=1; n<Ny; n++){
            int b = index+n;
            auto g = BondGate(sites,b,b+1);
            auto AA = psi_DM(b)*psi_DM(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_DM.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi_DM.position(g.i2()); //orthogonality center moves to the right
          }
        }//long-range interaction
      }// for j
    }// for i
    psi_DM.orthogonalize({"Cutoff=",1E-12}); //restore psi_DM bond dimensions
    
    //write to file
    enerfile1 << tval[i+1] << " " << energy_DM << " " << maxLinkDim(psi_DM) << " "; //print to file
    for(int j = 0; j<N; j++){ //save local energy values
      enerfile1 << LocalEnergy[j] << " ";
    }
    enerfile1 << std::endl;

    //check Fit method for MPO*MPS
    psi_Fit = applyMPO(expH1,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels
    psi_Fit = applyMPO(expH2,psi_Fit,args_Fit);
    psi_Fit.noPrime().normalize(); //need to do this after each to take care of prime levels
    auto energy_Fit = innerC(psi_Fit,H,psi_Fit).real();

    // calculate local energy density
    for(int i=1; i<=Nx; i++){
      for(int j=1; j<=Ny; j++){ 
        int index = (i-1)*Ny + j; //this order to make orthogonality center easier to compute
        ITensor ket;
        if(j==Ny){ //y-periodic boundary equations with swap gates
          psi_Fit.position(index-Ny+1);
          for(int n=1; n<Ny-1; n++){
            int b = index-Ny+n; //define gate bond
            auto g = BondGate(sites,b,b+1);
            auto AA = psi_Fit(b)*psi_Fit(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_Fit.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi_Fit.position(g.i2()); //orthogonality center moves to the right
          }
          
          ket = psi_Fit(index-1)*psi_Fit(index);
          LocalEnergy[index-1] = eltC( dag(prime(ket,"Site")) * LEDyPBC[i-1] * ket).real();

          //restore the state to the original MPS
          for(int n=Ny-1; n>1; n--){
            int b = index-Ny+n;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi_Fit(b-1)*psi_Fit(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0");
            psi_Fit.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi_Fit.position(g.i1()); //move orthogonality center to the left  
          }
        }// y-periodic

        //original nearest-neighbour code
        else{ 
          psi_Fit.position(index);
          ket = psi_Fit(index)*psi_Fit(index+1);
          LocalEnergy[index-1] = eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();
        } //nearest-neighbour

        if(i<Nx){
          psi_Fit.position(index+Ny); //site to bring to index+1  

          // bring index+Ny to position index+1
          for(int n=Ny; n>1; n--){
            int b = index+n;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi_Fit(b-1)*psi_Fit(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_Fit.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi_Fit.position(g.i1()); //orthogonality center moves to the left
          }

          ket = psi_Fit(index)*psi_Fit(index+1);
          LocalEnergy[index-1] += eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();

          // bring index+1 back to position index+Ny
          for(int n=1; n<Ny; n++){
            int b = index+n;
            auto g = BondGate(sites,b,b+1);
            auto AA = psi_Fit(b)*psi_Fit(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi_Fit.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi_Fit.position(g.i2()); //orthogonality center moves to the right
          }
        }//long-range interaction
      }// for j
    }// for i
    psi_Fit.orthogonalize({"Cutoff=",1E-12}); //restore psi_Fit bond dimensions
    
    //write to file
    enerfile2 << tval[i+1] << " " << energy_Fit << " " << maxLinkDim(psi_Fit) << " ";
    for(int j = 0; j<N; j++){ //save local energy values
      enerfile2 << LocalEnergy[j] << " ";
    }
    enerfile2 << std::endl;

    printfln("Iteration %d, DensityMatrix energy = %0.3g, max link dim is %d",i+1,energy_DM,maxLinkDim(psi_DM));
    printfln("            ,          Fit energy = %0.3g, max link dim is %d",energy_Fit,maxLinkDim(psi_Fit));

  }
  
  enerfile1.close(); enerfile2.close();

  return 0;
  }
