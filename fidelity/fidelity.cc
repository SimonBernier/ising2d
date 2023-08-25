//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[]){
    
    std::clock_t tStart = std::clock();

    if(argc < 2){
        printfln("Usage: %s input_file",argv[0]);
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto Ly = input.getInt("Ly", 3);
    auto Lx = input.getInt("Lx", 16);
    auto dh = input.getReal("dh", 0.01);
    auto hstep = input.getReal("hstep",0.01);
    auto h0 = input.getReal("h0",1.0);
    auto hmax = input.getReal("hmax",4.0);
    int ISTOP = int((hmax-h0)/hstep);
    
    printfln("Ly = %d, Lx = %d, dh = %0.4f", Ly, Lx, dh);

    //write results to file
    char schar[64];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_dh_%0.4f_2dTFI_fidelity.dat",Ly,Lx,dh);
    std::string s1(schar);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "h" << " " << "en0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "var1" << " " << "maxBondDim1" << " " << "1-F" <<  std::endl;


    auto N = Lx * Ly;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});
    
    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*h0, "Sz", j;
    }
    auto H = toMPO(ampo);

    //initial state
    auto state = InitState(sites); 
    for(auto j : range1(N)){
        state.set(j, "Up");
    }
    auto initState = MPS(state);

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 10, 20, 100, 100, 200, 200, 400, 400, 512;
    sweeps.cutoff() = 1E-10;

    //
    //solve for ground state of h-dh/2
    //
    // first h-dh state
    auto [en0,psi0] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = sqrt( abs( inner(H,psi0,H,psi0)-en0*en0) );
    auto maxBondDim = maxLinkDim(psi0);

    printfln("\nh = %0.2f", h0);
    println("h state");
    printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en0, var, maxBondDim);

    dataFile << h0 << " " << en0 << " " << var << " " << maxBondDim << " ";

    //
    // ground state of h+dh/2
    //
    for(auto j : range1(N)){
        ampo += -2.0*dh, "Sz", j;
    }
    H = toMPO(ampo);

    auto [en1,psi1] = dmrg(H,psi0,sweeps,{"Silent=",true});
    var = sqrt( abs( inner(H,psi1,H,psi1)-en1*en1) );
    maxBondDim = maxLinkDim(psi1);

    println("h+dh/2 state");
    printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en1, var, maxBondDim);

    dataFile << en1 << " " << var << " " << maxBondDim << " ";

    //
    // calculate fidelity
    //
    auto F = abs(inner(psi1,psi0));
    printfln("Fidelity F = %0.10f", F);

    dataFile << 1.-F << " " << std::endl;
    
    for(int i = 1; i <= ISTOP ; i++){
        double h = h0 + i*hstep;
        printfln("\nh = %0.2f", h);

        // autompo hamiltonian
        ampo = AutoMPO(sites);
        for(auto j : lattice){
            ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
        }
        for(auto j : range1(N)){
            ampo += -2.0*h, "Sz", j;
        }
        H = toMPO(ampo);

        en0 = dmrg(psi0,H,sweeps,{"Silent=",true});
        var = sqrt( abs( inner(H,psi0,H,psi0)-en0*en0) );
        maxBondDim = maxLinkDim(psi0);

        println("h state");
        printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en0, var, maxBondDim);

        dataFile << h << " " << en0 << " " << var << " " << maxBondDim << " ";

        //
        // ground state of h+dh
        //
        for(auto j : range1(N)){
            ampo += -2.0*dh, "Sz", j;
        }
        H = toMPO(ampo);

        en1 = dmrg(psi1,H,sweeps,{"Silent=",true});
        var = sqrt( abs( inner(H,psi1,H,psi1)-en1*en1) );
        maxBondDim = maxLinkDim(psi1);

        println("h+dh state");
        printfln("Energy = %0.3f, var = %0.3g, maxLinkDim = %d", en1, var, maxBondDim);

        dataFile << en1 << " " << var << " " << maxBondDim << " ";

        //
        // calculate fidelity
        //
        F = abs(inner(psi1,psi0));
        printfln("Fidelity F = %0.10f", F);

        dataFile << 1.-F << " " << std::endl;
    }

    dataFile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

}//main
