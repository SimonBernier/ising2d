#include "itensor/all.h"

using namespace itensor;

//function definition for calculation of local energy
std::vector<double> localEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);

int main(int argc, char *argv[])
  {
  int Ly = 4;
  int Lx = 16;
  float h = 4.0;

  if(argc > 3)
    h = std::stof(argv[3]);
  if(argc > 2)
    Lx = std::stoi(argv[2]);
  if(argc > 1)
    Ly = std::stoi(argv[1]);

  // write results to file
  char schar[64];
  int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.3g_CP.dat",Ly,Lx,h); 
  std::string s1(schar);
  std::ofstream enerfile;
  enerfile.open(s1); // opens the file
  if( !enerfile ) { // file couldn't be opened
      std::cerr << "Error: file could not be opened" << std::endl;
      exit(1);
  }
  //make header
  enerfile << "energy" << " " << "var" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;

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
  sweeps.maxdim() = 20, 50, 100, 200, 400, 800;
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
        LEDyPBC[i-1] += -2.0*h*sites.op("Id",index-1)*sites.op("Sx",index);
      }
      // MPS nearest-neighbour
      if(j<Ly){
        LED[i-1][j-1] = -4.0*sites.op("Sz",index)*sites.op("Sz",index+1);
        LED[i-1][j-1] += -2.0*h*sites.op("Sx",index)*sites.op("Id",index+1);
      }
    }
  }

  // calculate ground state of critical H
  auto [energy, psi] = dmrg(H,MPS(state),sweeps,{"Silent=",true});
  auto var = inner(psi,H,H,psi) - energy*energy;
  // calculate local energy
  LocalEnergy = localEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);
  // store to file
  enerfile << energy << " " << var << " " << maxLinkDim(psi) << " ";
  for(int j = 0; j<L; j++){ //save local energy values
  enerfile << LocalEnergy[j] << " ";
  }
  enerfile << std::endl;

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
