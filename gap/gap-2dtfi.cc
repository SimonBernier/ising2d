#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
    {
    
    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);
    
    int Ly, Lx;
    double h;

    std::ifstream parameter_file ("parameter_list.txt");
    std::string parameter;
    int temp = 0;
    if ( parameter_file.is_open() ) { // always check whether the file is open
        while(temp<=runNumber){ //skip the appropriate number of lines
            std::getline(parameter_file, parameter);
            temp += 1;
        } 
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        Ly = std::stoi(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        Lx = std::stoi(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        h = std::stod(parameter);
    }
    parameter_file.close();
    
    printfln("Ly = %d, Lx = %d, h = %0.3f", Ly, Lx, h);

    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.3f_2dTFI_gap.dat",Ly,Lx,h);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "en0" << " " << "M0" << " " << "p0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "M1" << " " << "p1" << " " << "var1" << " " << "maxBondDim1" << " "
             << "en2" << " " << "M2" << " " << "p2" << " " << "var2" << " " << "maxBondDim2" << " " << std::endl;

    auto N = Lx * Ly;
    auto sites = SpinHalf(N,{"ConserveSz=",false});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});
    
    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sz", j.s1, "Sz", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*h, "Sx", j;
    }

    // make parity operator
    auto P = MPO(N);
    for(auto j : range1(N)){
        P.ref(j) = 2.0*sites.op("Sx",j);
    }

    // make vector of Sz operators
    std::vector<ITensor> Sz(N);
    for(auto j : range1(N)){
        Sz[j-1] = 2.0*sites.op("Sz",j);
    }

    // 2d ising model parameters
    auto sweeps = Sweeps(15);
    sweeps.maxdim() = 10, 20, 100, 100, 200, 200, 400, 400, 512;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 1E-7,1E-8,0.0;

    //
    //solve for ground state
    //
    auto H = toMPO(ampo);
    auto [en0,psi0] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
    auto var = inner(psi0,H,H,psi0)-en0*en0;
    auto maxBondDim = maxLinkDim(psi0);
    double M0 = 0.0;
    for(auto b : range1(N)){
        psi0.position(b);
        auto m = elt( dag(prime(psi0(b),"Site")) * Sz[b-1] * psi0(b) );
        M0 += m;
    }
    M0 /= double(N);
    auto p0 = inner(psi0,P,psi0);
    println("\nfirst state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en0, M0, p0, maxBondDim, var);

    dataFile << en0 << " " << M0 << " " << p0 << " " << var << " " << maxBondDim << " ";

    //
    // excited state for PM regime, second GS for FM regime
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    var = inner(psi1,H,H,psi1)-en1*en1;
    maxBondDim = maxLinkDim(psi1);
    double M1 = 0.0;
    for(auto b : range1(N)){
        psi1.position(b);
        auto m = elt( dag(prime(psi1(b),"Site")) * Sz[b-1] * psi1(b) );
        M1 += m;
    }
    M1 /= double(N);
    auto p1 = inner(psi1,P,psi1);
    println("\nsecond state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en1, M1, p1, maxBondDim, var);

    dataFile << en1 << " " << M1 << " " << p1 << " " << var << " " << maxBondDim << " ";

    //
    // first excited state in the ferromagnetic regime
    //
    wfs.push_back(psi1);
    auto [en2,psi2] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    var = inner(psi2,H,H,psi2)-en2*en2;
    maxBondDim = maxLinkDim(psi2);
    double M2 = 0.0;
    for(auto b : range1(N)){
        psi2.position(b);
        auto m = elt( dag(prime(psi2(b),"Site")) * Sz[b-1] * psi2(b) );
        M2 += m;
    }
    M2 /= double(N);
    auto p2 = inner(psi2,P,psi2);
    println("\nthird state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en2, M2, p2, maxBondDim, var);

    dataFile << en2 << " " << M2 << " " << p2 << " " << var << " " << maxBondDim << " " << std::endl;


    dataFile.close();

    return 0;
}// runPT

