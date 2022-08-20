#include "itensor/all.h"

using namespace itensor;

//function definition for calculation of local energy
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);

double calculatePBCenergy(int i, int Ly, MPS psi, ITensor LED, SiteSet sites);

std::vector<double> calculateLRenergy(int i, int Ly, MPS psi, std::vector<std::vector<ITensor>> LED_LR, SiteSet sites);

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[]){

    int Ly = 2;
    int Lx = 16;
    float h = 4.0;

    if(argc > 3)
        h = std::stof(argv[3]);
    if(argc > 2)
        Lx = std::stoi(argv[2]);
    if(argc > 1)
        Ly = std::stoi(argv[1]);

    printfln("Ly = %d, Lx = %d, h = %0.2f", Ly, Lx, h);

    // write results to file
    char schar[64];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_tfi2Dcrit_En.dat",Ly,Lx,h); 
    std::string s1(schar);
    std::ofstream enerfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
          std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    enerfile << "energy" << " " << "var" << " " << "SvN" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;

    auto N = Ly*Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*h, "Sz", j;
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
    auto sweeps = Sweeps(15);
    sweeps.maxdim() = 10, 20, 100, 100, 200, 200, 400, 400, 512;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 1E-7,1E-8,0.0;

    // calculate initial local energy density
    std::vector<double> localEnergy(Ly*(Lx-1),0.0); // local energy density vector
  
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
                LED_LR[i-1][m] = -4.0*sites.op("Sx",index+2*m)*sites.op("Sx",index+2*m+1);
                }
            }
            //y-periodic boundary equations
            if(j==Ly){
                // site index-Ly+1 is moved to site index-1 with swap gates
                LEDyPBC[i-1] = -4.0*sites.op("Sx",index-1)*sites.op("Sx",index);
                LEDyPBC[i-1] += -h *sites.op("Id",index-1)*sites.op("Sz",index);
                LEDyPBC[i-1] += -h *sites.op("Sz",index-1)*sites.op("Id",index);
            }
            // MPS nearest-neighbour
            if(j<Ly){
                LED[i-1][j-1] = -4.0*sites.op("Sx",index)*sites.op("Sx",index+1);
                LED[i-1][j-1] += -h *sites.op("Sz",index)*sites.op("Id",index+1);
                LED[i-1][j-1] += -h *sites.op("Id",index)*sites.op("Sz",index+1);
            }
        }
    }

    // calculate ground state of critical H
    auto [energy, psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H,psi,H,psi) - energy*energy;
    // calculate local energy
    localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);
    // calculate von Neumann entanglement entropy
    auto svN = vonNeumannS(psi, N/2);

    printfln("\n energy = %0.2f, var = %0.2g, SvN = %0.3f, bondDim = %d", energy, var, svN, maxLinkDim(psi));

    // store to file
    enerfile << energy << " " << var << " " << svN << " " << maxLinkDim(psi) << " ";
    for(int j = 0; j<Ly*(Lx-1); j++){ //save local energy values
    enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;

    enerfile.close();

    return 0;
    }

// calculates local energy for 2D MPS using gates
std::vector<double> calculateLocalEnergy(int Lx, int Ly, SiteSet sites, MPS psi,
                                std::vector<std::vector<ITensor>> LED, 
                                std::vector<ITensor> LEDyPBC,
                                std::vector<std::vector<ITensor>> LED_LR){
  
    std::vector<std::vector<double>> tempEn(Lx, std::vector<double>(Ly, 0.0));

    std::vector<double> localEnergy(Ly*(Lx-1), 0.0); // interpolated energy density

    // MPS nearest-neighbour interaction
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<Ly; j++){ 

            int index = (i-1)*Ly + j;

            psi.position(index);
            auto ket = psi(index)*psi(index+1);
            tempEn[i-1][j-1] = eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();

        }// for j
    }// for i

    // y-periodic interactions
    for(int i=1; i<=Lx; i++){

        tempEn[i-1][Ly-1] = calculatePBCenergy(i, Ly, psi, LEDyPBC[i-1], sites);

    }// for i

    // MPS long-range interactions
    for(int i=1; i<Lx; i++){

        int index = (i-1)*Ly + 1;
        auto lrEnergy = calculateLRenergy(i, Ly, psi, LED_LR, sites);

        for(int m = 0; m<Ly; m++){

            localEnergy[index-1+m] = lrEnergy[m];

        } // for m

    }// for i

    // interpolate tempEn into localEnergy
    for(int i=1; i<Lx; i++){
        for(int j=1; j<=Ly; j++){ 

            int index = (i-1)*Ly + j;

            if(j == 1){ // interpolate with periodic boundary conditions
                localEnergy[index-1] = 0.25 * ( tempEn[i-1][j-1] + tempEn[i][j-1] + tempEn[i-1][j+Ly-2] + tempEn[i][j+Ly-2] );
            }
            else{ // interpolate normally
                localEnergy[index-1] = 0.25 * ( tempEn[i-1][j-2] + tempEn[i][j-2] + tempEn[i-1][j-1] + tempEn[i][j-1] );
            }

        } // for j
    } // for i

    return localEnergy;

}//localEnergy

double calculatePBCenergy(int i, int Ly, MPS psi, ITensor hterm, SiteSet sites){

    double energy;
    int index = (i-1)*Ly + 1;

    psi.position(index);

    for(int n=0; n<Ly-2; n++){

        int b = index+n;//define gate bond
        auto g = BondGate(sites,b,b+1);
        auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
        AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
        psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
        psi.position(g.i2()); //orthogonality center moves to the right

        //printf("sg (%d,%d)", b,b+1);

    } // for n

    //println("");
                
    auto ket = psi(index+Ly-2)*psi(index+Ly-1);
    energy = eltC( dag(prime(ket,"Site")) * hterm * ket).real();

    //printfln("e(%d,%d)", index+Ly-2, index+Ly-1);

    /*
    //restore the state to the original MPS
    for(int n=Ly-2; n>0; n--){
        int b = index+n;
        auto g = BondGate(sites,b-1,b);
        auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
        AA.replaceTags("Site,1","Site,0");
        psi.svdBond(g.i1(), AA, Fromright); //svd from the right
        psi.position(g.i1()); //move orthogonality center to the left  
    } // for n
    */
    
    return energy;
} // calculatePBCenergy

std::vector<double> calculateLRenergy(int i, int Ly, MPS psi, std::vector<std::vector<ITensor>> LED_LR, SiteSet sites){

    std::vector<double> energy(Ly,0.);

    int index = (i-1)*Ly + 1;

    //smart ordering of gates 
    for(int m=0; m<=Ly-2; m++){

        psi.position(index+Ly+m);

        for(int n=Ly; n>m+1; n--){

            int b = index + n + m;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi.position(g.i1()); //orthogonality center moves to the left

            //printf("sg (%d,%d)", b-1,b);

        } // for n

        //println("");

    } // for m
                
    for(int m = 0; m<Ly; m++){

        psi.position(index+2*m);
        auto ket = psi(index+2*m)*psi(index+2*m+1);
        energy[m] = eltC( dag(prime(ket,"Site")) * LED_LR[i-1][m] * ket).real();

        //printf("e(%d,%d)", index+2*m, index+2*m+1);
    } // for m
    //println("");
    /*
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
        } // for n
    }// for m
    */

    return energy;
}

//calculate entanglement
Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //choose orthogonality center and perform svd
    psi.position(b);
    auto l = leftLinkIndex(psi,b);
    auto s = siteIndex(psi,b);
    auto [U,S,V] = svd(psi(b),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    for(auto n : range1(dim(u))){
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    return SvN;

}//vonNeumannS