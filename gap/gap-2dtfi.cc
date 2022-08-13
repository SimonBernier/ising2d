#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int Lx=16, Ly=3;
    double h=1.;

    if(argc > 2)
        Lx = std::stoi(argv[2]);
    if(argc > 3)
        Ly = std::stoi(argv[3]);
    if(argc > 4)
        h = std::stod(argv[4]);   
    
    printfln("Ly = %d, Lx = %d, h = %0.2f", Ly, Lx, h);

    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.2f_2dTFI_gap.dat",Ly,Lx,h);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "en0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "var1" << " " << "maxBondDim1" << " " << std::endl;

    auto N = Lx * Ly;
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

    //
    //solve for ground state
    //
    auto [en0,psi0] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H,psi0,H,psi0) - en0*en0;
    auto maxBondDim = maxLinkDim(psi0);
    println("\nfirst state");
    printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en0, var, maxBondDim);

    dataFile << en0 << " " << var << " " << maxBondDim << " ";

    //
    // excited state
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(H,wfs,initState,sweeps,{"Silent=",true,"Weight=",20.0});
    var = inner(H,psi1,H,psi1) - en1*en1;
    maxBondDim = maxLinkDim(psi1);

    println("\nsecond state");
    printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en1, var, maxBondDim);

    dataFile << en1 << " " << var << " " << maxBondDim << " ";

    dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}// runPT

