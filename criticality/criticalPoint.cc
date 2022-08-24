#include "itensor/all.h"

using namespace itensor;

//std::vector<MPO> makeEnergyMPO(int, int, double, MPS, SiteSet);

//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);
// calculate energy for periodic boundary at x = i
double calculatePBCenergy(int, int, MPS, ITensor, SiteSet);
// calculate energy of horizontal bonds for one column
std::vector<double> calculateLRenergy(int, int, MPS, std::vector<std::vector<ITensor>>, SiteSet);

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
double spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

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
    sweeps.noise() = 1E-7,1E-8,1E-7,1E-8,1E-7,1E-8,1E-7,1E-8,1E-7,1E-8,0;

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
                    LED_LR[i-1][m] += -h *sites.op("Sz",index+2*m)*sites.op("Id",index+2*m+1);
                    LED_LR[i-1][m] += -h *sites.op("Id",index+2*m)*sites.op("Sz",index+2*m+1);
                    if( i==1 ){ // add to left
                        LED_LR[i-1][m] += -h*sites.op("Sz",index+2*m)*sites.op("Id",index+2*m+1);
                    }
                    else if( i==Lx-1){ // add to right
                        LED_LR[i-1][m] += -h*sites.op("Id",index+2*m)*sites.op("Sz",index+2*m+1);
                    }
                }
            }
            //y-periodic boundary equations
            if(j==Ly){
                // site index-Ly+1 is moved to site index-1 with swap gates
                LEDyPBC[i-1] = -4.0*sites.op("Sx",index-1)*sites.op("Sx",index);
            }
            // MPS nearest-neighbour
            if(j<Ly){
                LED[i-1][j-1] = -4.0*sites.op("Sx",index)*sites.op("Sx",index+1);
            }
        }
    }

    // Create the Local Energy Density MPOs
    //std::vector<MPO> LEDMPO = makeEnergyMPO(Lx,Ly,h,psi,sites);

    // calculate ground state of critical H
    auto [energy, psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H,psi,H,psi) - energy*energy;

    // calculate local energy
    //for( int i = 1; i <= (Lx-1)*Ly; i++){
    //    localEnergy[i-1] = inner(psi,LED[i-1],psi);
    //}
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

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    }

/*
std::vector<MPO> makeEnergyMPO(int Lx, int Ly, double h, MPS psi, SiteSet sites){
    
    std::vector<MPO> LED( (Lx-1)*Ly );
    for (int i = 1; i < Lx; i++){
        for (int j = 1; j<=Ly; j++){

            int b = (i-1)*Ly + j;
            auto ampo = AutoMPO(sites);

            // long-range
            ampo += -4.0, "Sx", b, "Sx", b+Ly;

            // on-site transverse field
            if(i == 1){ // add missing terms at the left side of the lattice
                ampo += -2.0*h , "Sz", b;
                ampo += -h , "Sz", b+Ly;
                if(j==1){
                    ampo += -1.0, "Sx", b, "Sx", b+1;
                    ampo += -1.0, "Sx", b, "Sx", b+Ly-1;
                }
                else if(j == Ly){
                    ampo += -1.0, "Sx", b-1, "Sx", b;
                    ampo += -1.0, "Sx", b-Ly+1, "Sx", b;
                }
                else{
                    ampo += -1.0, "Sx", b, "Sx", b+1;
                    ampo += -1.0, "Sx", b-1, "Sx", b;
                }
            }
            else if(i == Lx-1){ // add missing terms at the right side of the lattice
                ampo += -h , "Sz", b;
                ampo += -2.0*h , "Sz", b+Ly;
                if(j==1){
                    ampo += -1.0, "Sx", b+Ly, "Sx", b+Ly+1;
                    ampo += -1.0, "Sx", b+Ly, "Sx", b+2*Ly-1;
                }
                else if(j == Ly){
                    ampo += -1.0, "Sx", b+Ly-1, "Sx", b+Ly;
                    ampo += -1.0, "Sx", b+1, "Sx", b+Ly;
                }
                else{
                    ampo += -1.0, "Sx", b+Ly, "Sx", b+Ly+1;
                    ampo += -1.0, "Sx", b+Ly-1, "Sx", b+Ly;
                }
            }
            else{ // center of the chain
                ampo += -h , "Sz", b;
                ampo += -h , "Sz", b+Ly;
            }
            
            // interpolation
            if( j == Ly ){ //top of the lattice
                ampo += -1.0, "Sx", b-1, "Sx", b;
                ampo += -1.0, "Sx", b, "Sx", b-Ly+1;
                ampo += -1.0, "Sx", b+Ly-1, "Sx", b+Ly;
                ampo += -1.0, "Sx", b+Ly, "Sx", b+1;
            }
            else if( j == 1){ //bottom of the lattice
                ampo += -1.0, "Sx", b+Ly-1, "Sx", b;
                ampo += -1.0, "Sx", b, "Sx", b+1;
                ampo += -1.0, "Sx", b+2*Ly-1, "Sx", b+Ly;
                ampo += -1.0, "Sx", b+Ly, "Sx", b+Ly+1;
            }
            else{ // middle of the lattice
                ampo += -1.0, "Sx", b-1, "Sx", b;
                ampo += -1.0, "Sx", b, "Sx", b+1;
                ampo += -1.0, "Sx", b+Ly-1, "Sx", b+Ly;
                ampo += -1.0, "Sx", b+Ly, "Sx", b+Ly+1;
            }
            
            // make the MPO
            LED[b-1] = toMPO(ampo);
            //printfln("max link dimension of MPO at bond %d = %d", b, maxLinkDim(LED[b-1]));
        }// for j
    }// for i

    return LED;
}
*/

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
                localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-1] + tempEn[i][j-1] + tempEn[i-1][j+Ly-2] + tempEn[i][j+Ly-2] );
            }
            else{ // interpolate normally
                localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-2] + tempEn[i][j-2] + tempEn[i-1][j-1] + tempEn[i][j-1] );
            }
            // add left/right boundary terms
            if( i==1 ){
                if( j== 1){
                    localEnergy[index-1] += 0.25 * ( tempEn[i-1][j+Ly-2] + tempEn[i-1][j-1]);
                }
                else{
                    localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-2] + tempEn[i-1][j-1]);
                }
            }// if i=1
            else if( i==Lx-1 ){
                if( j== 1){
                    localEnergy[index-1] += 0.25 * ( tempEn[i][j+Ly-2] + tempEn[i][j-1]);
                }
                else{
                    localEnergy[index-1] += 0.25 * ( tempEn[i][j-2] + tempEn[i][j-1]);
                }
            } // if i=Lx-1
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

    //restore the state to the original MPS
    for(int n=Ly-2; n>0; n--){
        int b = index+n;
        auto g = BondGate(sites,b-1,b);
        auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
        AA.replaceTags("Site,1","Site,0");
        psi.svdBond(g.i1(), AA, Fromright); //svd from the right
        psi.position(g.i1()); //move orthogonality center to the left

        //printf("sg (%d,%d) ", b-1,b);

    } // for n
    
    //println("\n");
    
    return energy;
} // calculatePBCenergy

std::vector<double> calculateLRenergy(int i, int Ly, MPS psi, std::vector<std::vector<ITensor>> LED_LR, SiteSet sites){

    std::vector<double> energy(Ly,0.);

    int index = (i-1)*Ly + 1;

    // long range interactions with smart ordering of gates 
    for(int m=0; m<=Ly-2; m++){

        psi.position(index+Ly+m);

        for(int n=Ly; n>m+1; n--){

            int b = index + n + m;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi.position(g.i1()); //orthogonality center moves to the left

            //printf("sg (%d,%d) ", b-1,b);

        } // for n

        //println("");

    } // for m
                
    for(int m = 0; m<Ly; m++){

        psi.position(index+2*m);
        auto ket = psi(index+2*m)*psi(index+2*m+1);
        energy[m] = eltC( dag(prime(ket,"Site")) * LED_LR[i-1][m] * ket).real();

        //printf("e(%d,%d) ", index+2*m, index+2*m+1);
    } // for m
    
    //println("");

    // bring index+1 back to position index+Ly
    for(int m=Ly-2; m>=0; m--){

        psi.position(index+1+2*m);

        for(int n=m+1; n<Ly; n++){

            int b = index + n + m ;
            auto g = BondGate(sites,b,b+1);
            auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
            psi.position(g.i2()); //orthogonality center moves to the right

            //printf("sg (%d,%d) ", b,b+1);

        } // for n

        //println("");

    }// for m

    //println("");
    
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

//calculate spin-spin correlator
double spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrX;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left 
        }
        auto ket = psi(center)*psi(center+1);
        auto SxSx = 4.0*sites.op("Sx",center)*sites.op("Sx",center+1);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SxSx = 4.0*sites.op("Sx",center-1)*sites.op("Sx",center);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else{
        corrX = 1.;
    }

    return corrX;

}//SxSx
