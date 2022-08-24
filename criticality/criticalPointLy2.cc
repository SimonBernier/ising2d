#include "itensor/all.h"

using namespace itensor;

//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS, std::vector<ITensor>,
                                                                 std::vector<std::vector<ITensor>>);
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
    sweeps.noise() = 1E-7,1E-8,0;

    // calculate initial local energy density
    std::vector<double> localEnergy((Lx-1)*Ly,0.0); // local energy density vector
  
    //make 2D vector of ITensor for local energy operators
    //long-range interactions have the same structure as nearest-neighbour when we use swap gates
    std::vector<ITensor> LED(Lx);
    std::vector<std::vector<ITensor>> LED_LR(Lx-1, std::vector<ITensor>(Ly));
    // make local energy tensors
    for(int i=1; i<=Lx; i++){

        int b = (i-1)*Ly+1;

        // MPS nearest-neighbour
        LED[i-1] = -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);

        //MPS long-range
        if(i<Lx){
            LED_LR[i-1][0] = -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
            LED_LR[i-1][1] = -4.0*sites.op("Sx",b+Ly)*sites.op("Sx",b+Ly+1);
            LED_LR[i-1][0] += -h*sites.op("Sz",b)*sites.op("Id",b+1);
            LED_LR[i-1][0] += -h*sites.op("Id",b)*sites.op("Sz",b+1);
            LED_LR[i-1][1] += -h*sites.op("Sz",b+Ly)*sites.op("Id",b+Ly+1);
            LED_LR[i-1][1] += -h*sites.op("Id",b+Ly)*sites.op("Sz",b+Ly+1);
            if( i==1 ){ // add to left
                LED_LR[i-1][0] += -h*sites.op("Sz",b)*sites.op("Id",b+1);
                LED_LR[i-1][1] += -h*sites.op("Sz",b+Ly)*sites.op("Id",b+Ly+1);
            }
            else if( i==Lx-1){ // add to right
                LED_LR[i-1][0] += -h*sites.op("Id",b)*sites.op("Sz",b+1);
                LED_LR[i-1][1] += -h*sites.op("Id",b+Ly)*sites.op("Sz",b+Ly+1);
            }
        }
    }

    // calculate ground state of critical H
    auto [energy, psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H,psi,H,psi) - energy*energy;

    // calculate local energy
    localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LED_LR);
    // calculate von Neumann entanglement entropy
    auto svN = vonNeumannS(psi, N/2);

    printfln("\n energy = %0.2f, var = %0.2g, SvN = %0.3f, bondDim = %d", energy, var, svN, maxLinkDim(psi));

    // store to file
    enerfile << energy << " " << var << " " << svN << " " << maxLinkDim(psi) << " ";
    for(int j = 0; j<(Lx-1)*Ly; j++){ //save local energy values
    enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;

    enerfile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    }

// calculates local energy for 2D MPS using gates
std::vector<double> calculateLocalEnergy(int Lx, int Ly, SiteSet sites, MPS psi,
                                std::vector<ITensor> LED, 
                                std::vector<std::vector<ITensor>> LED_LR){
  
    std::vector<double> vertBonds(Lx, 0.0);
    std::vector<std::vector<double>> horBonds((Lx-1), std::vector<double>(Ly, 0.0)); 
    std::vector<double> localEnergy( (Lx-1)*Ly, 0.0); // interpolated energy density

    // MPS nearest-neighbour interaction
    for(int i=1; i<=Lx; i++){
        int index = (i-1)*Ly + 1;

        psi.position(index);
        auto ket = psi(index)*psi(index+1);
        vertBonds[i-1] = eltC(dag(prime(ket,"Site")) * LED[i-1] * ket).real();

    }// for i

    // MPS long-range interactions
    for(int i=1; i<Lx; i++){

        auto lrEnergy = calculateLRenergy(i, Ly, psi, LED_LR, sites);

        for(int m = 0; m<Ly; m++){

            horBonds[i-1][0] = lrEnergy[0];
            horBonds[i-1][1] = lrEnergy[1];

        } // for m

    }// for i
    
    // interpolate into localEnergy
    for(int i=1; i<Lx; i++){
        for(int j=1; j<=Ly; j++){
            int b = (i-1)*Ly + j;
            localEnergy[b-1] = horBonds[i-1][j-1];
            localEnergy[b-1] += 0.25 * ( vertBonds[i-1] + vertBonds[i] );
            if( i==1 ){
                localEnergy[b-1] += 0.25 * vertBonds[i-1];
            }
            else if( i==Lx-1 ){
                localEnergy[b-1] += 0.25 * vertBonds[i];
            }
        }
    } // for i

    return localEnergy;

}//localEnergy

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

        } // for n

    } // for m
                
    for(int m = 0; m<Ly; m++){

        psi.position(index+2*m);
        auto ket = psi(index+2*m)*psi(index+2*m+1);
        energy[m] = eltC( dag(prime(ket,"Site")) * LED_LR[i-1][m] * ket).real();

    } // for m
    
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

        } // for n

    }// for m
    
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
