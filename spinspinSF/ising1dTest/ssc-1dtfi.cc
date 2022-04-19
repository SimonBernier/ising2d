#include "itensor/all.h"

using namespace itensor;

//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, double, SiteSet, std::vector<ITensor>);

int main(int argc, char *argv[])
  {
  int N = 16;
  float h = 4.0;
  int Nt = 200;
  float dt = 0.1;

  if(argc > 2)
    h = std::stof(argv[2]);
  if(argc > 1)
    N = std::stoi(argv[1]);

  printfln("N = %d, h = %0.2f", N,h);

  // write results to file
  char schar[64];
  int n1 = std::sprintf(schar,"N_%d_h_%0.3g_sscTFI1d.dat",N,h);

  std::string s1(schar);
  std::ofstream dataFile;
  dataFile.open(s1); // opens the file
  if( !dataFile ) { // file couldn't be opened
    std::cerr << "Error: file could not be opened" << std::endl;
    exit(1);
  }
  //make header for t=0 calculations
  dataFile << "t=0" << " " << "enPsi" << " " << "MaxDimPsi" << " " << "enPhi" << " " << "MaxDimPhi" << " " << "Sz(x)Sz(N/2,0)" << " " << std::endl;

  auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

  auto ampo = AutoMPO(sites);

  // autompo hamiltonian
  for (int b = 1; b < N; b++){
    ampo += -4.0,"Sx", b, "Sx", b+1;
  }
  for (int b = 1; b <= N; b++){
    ampo += -2.0,"Sz",b;
  }
  auto H = toMPO(ampo);

  //initial state
  auto state = InitState(sites); 
  for(int i = 1; i <= N; i++){
    state.set(i,"Up");
  }
  auto initState = MPS(state);
  PrintData(totalQN(initState));

  // 2d ising model parameters
  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 20, 50, 100, 200;
  sweeps.cutoff() = 1E-12;
  
  //make vector of ITensor for local energy operators
  std::vector<ITensor> LED(N);
  // make local energy tensors
  for (int b = 1; b < N; b++){
    LED[b-1] =  -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
    LED[b-1] += -2.0*sites.op("Sz",b)*sites.op("Id",b+1);
  }
  LED[N-1] = -2.0*sites.op("Id",N-1)*sites.op("Sz",N);

  //make vector of Sz MPOs
  std::vector<MPO> Sz(N);
  for(auto j : range1(N)){
    auto ampoSz = AutoMPO(sites);
    for(auto i : range1(N)){
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
  auto [en_psi, psi] = dmrg(H,initState,sweeps,{"Silent=",true});

  // make |phi> = Sz|psi>
  int loc = N/2; //centered in x
  psi.position(loc);
  auto newA = 2.0*sites.op("Sz",loc) * psi(loc);
  newA.noPrime();
  auto phi = psi;
  phi.set(loc, newA);
  auto en_phi = inner(phi,H,phi);

  // store spin-spin correlation function
  std::vector<Complex> szsz(N,0.0); // local energy density vector
  //calculate the spin-spin correlation using MPS * MPO * MPS methods
  for(auto j : range1(N)){
    szsz[j-1] = innerC(psi, Sz[j-1], phi);
  }

  printfln("\nIteration %d, time = %0.2f; phi energy = %0.f, max link dim is %d",0,0, en_phi,maxLinkDim(phi));
  // store to file
  dataFile << 0.0 << " " << en_psi << " " << maxLinkDim(psi) << " " << en_phi << " " << maxLinkDim(phi) << " ";
  for(int j = 0; j<N; j++){ //save local energy values
    dataFile << real(szsz[j]) << " " << imag(szsz[j]) << " ";
  }
  dataFile << std::endl;

  dataFile << "tval" << " " << "enPhi" << " " << "MaxDimPhi" << " " << "Sz(x,t)Sz(N/2,0)" << " " << std::endl;

  // time evolution parameters
  double tval = 0.0; //time
  Real delta1 =  0.414490771794376*dt; // 1/(4-4^1/3)*dt
  Real delta2 = -0.657963087177503*dt; // (1-4*delta1/dt)*dt

  //args for time evolution methods
  Args args = Args("Cutoff=",1E-12,"MaxDim=",512);

  //Create a std::vector (dynamically sizeable array) to hold the Trotter gates
  std::vector<BondGate> gates;
  auto gatesdelta1 = makeGates(N, delta1, sites, LED);
  auto gatesdelta2 = makeGates(N, delta2, sites, LED);
  gates = gatesdelta1;
  gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
  gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
  gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
  gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));

  ////////////////////////////////////////////////////////////////////////////
  ///////// time evolve //////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  for(int n=1; n<=Nt; n++){
    tval += dt; //update time vector
    
    // do a time step
    gateTEvol(gates,dt,dt,phi,{args,"Verbose=",false});

    en_phi = innerC(phi, H, phi).real();

    //calculate the spin-spin correlation using MPS * MPO * MPS methods
    for(auto j : range1(N)){
      szsz[j-1] = innerC(exp(1_i*en_psi*tval)*psi, Sz[j-1], phi);
    }

    //write to file
    dataFile << tval << " " << en_phi << " " << maxLinkDim(phi) << " ";
    for(int j = 0; j<N; j++){ //save local energy values
      dataFile << real(szsz[j]) << " " << imag(szsz[j]) << " ";
    }
    dataFile << std::endl;

    printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(phi));

  }

  dataFile.close();

  return 0;
  }

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int L, double dt, SiteSet sites, std::vector<ITensor> LED)
  {
  std::vector<BondGate> gates; 
  //Create the gates exp(-i*tstep/2*hterm)
  for(int j=1; j<=L; j++){

    if(j<L){ // h_term = -J*Sz(j)Sz(j+1) - h*Sx(j)
      auto hterm = LED[j-1];
      auto g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm);
      gates.push_back(g);
    }
    else{ // h_term = - h*Sx(L)
      auto hterm = LED[j-1];
      auto g = BondGate(sites,j-1,j,BondGate::tReal,dt/2.,hterm);
      gates.push_back(g);  
    }//nearest-neighbour

  }// for j

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int j=L; j>=1; j--){

    //original nearest-neighbour code
    if(j==L){ // h_term = - h*Sx(L)
      auto hterm = LED[j-1];
      auto g = BondGate(sites,j-1,j,BondGate::tReal,dt/2.,hterm);
      gates.push_back(g);
    }
    else{ // h_term = -J*Sz(j)Sz(j+1) - h*Sx(j)
      auto hterm = LED[j-1];
      auto g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm);
      gates.push_back(g);
    } //nearest-neighbour

  }// for j

  return gates;
  
}// makeGates