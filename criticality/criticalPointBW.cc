//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"

using namespace itensor;

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

    printfln("Ly = %d, Lx = %d, h = %0.4f", Ly, Lx, h);
    // 2.5463 2.7300 2.8247 2.8806 2.9168 2.9418
    // 2.5400 2.7249 2.8202 2.8765 2.9130 2.9381

    // write results to file
    char schar[64];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.4f_tfi2DcritBW.dat",Ly,Lx,h); 
    std::string s1(schar);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
          std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    dataFile << "energy" << " " << "var" << " " << std::endl;

    auto N = Ly*Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity=",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += 4.0, "Sx", j.s1, "Sx", j.s2;
    }
    for(auto j : range1(N)){
        ampo += 2.0*h, "Sz", j;
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

    // calculate ground state of critical H
    auto [energy, psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = sqrt( abs( inner(H,psi,H,psi) - energy*energy) );

    printfln("energy = %0.5f, var = %0.5f, maxDim = %d", -energy, var, maxLinkDim(psi));

    // store to file
    dataFile << -energy << " " << var << " " << std::endl;

    dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    }