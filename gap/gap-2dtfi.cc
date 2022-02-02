#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
    {
    
    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);
    
    int Ly, Lx;
    double h, lambda;

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
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        lambda = std::stod(parameter);
    }
    parameter_file.close();
    
    printfln("Ly = %d, Lx = %d, h = %0.3f, lambda = %0.1f", Ly, Lx, h, lambda);

    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.3f_lambda_%0.1f_2dTFI_gap.dat",Ly,Lx,h,lambda);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "en0" << " " << "p0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "p1" << " " << "var1" << " " << "maxBondDim1" << " " << std::endl;

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
    auto P = MPO(sites);
    for(auto j : range1(N)){
        if(j==1){
            auto Pj = P(j);
            Pj.set(1,1,1, 0.0);
            Pj.set(1,2,1, 1.0);
            Pj.set(2,1,1, 1.0);
            Pj.set(2,2,1, 0.0);
            P.set(j,Pj);
        }
        else if(1<j && j<N){
            auto Pj = P(j);
            Pj.set(1,1,1,1, 0.0);
            Pj.set(1,2,1,1, 1.0);
            Pj.set(2,1,1,1, 1.0);
            Pj.set(2,2,1,1, 0.0);
            P.set(j,Pj);
        }
        else if(j==N){
            auto Pj = P(j);
            Pj.set(1,1,1, 0.0);
            Pj.set(1,2,1, 1.0);
            Pj.set(2,1,1, 1.0);
            Pj.set(2,2,1, 0.0);
            P.set(j,Pj);
        }        
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
    auto HP = H.plusEq(-lambda*P);
    auto [en0,psi0] = dmrg(HP,randomMPS(sites),sweeps,{"Silent=",true});
    en0 = inner(psi0,H,psi0);
    auto var = inner(psi0,H,H,psi0)-en0*en0;
    auto maxBondDim = maxLinkDim(psi0);
    auto p0 = inner(psi0,P,psi0);
    println("\nfirst state");
    printfln("Energy = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en0, p0, maxBondDim, var);

    dataFile << en0 << " " << p0 << " " << var << " " << maxBondDim << " ";

    //
    // excited state for FM regime, second ES for PM regime
    //
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;
    auto [en1,psi1] = dmrg(HP,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    en1 = inner(psi1,H,psi1);
    var = inner(psi1,H,H,psi1)-en1*en1;
    maxBondDim = maxLinkDim(psi1);
    auto p1 = inner(psi1,P,psi1);
    println("\nsecond state");
    printfln("Energy = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en1, p1, maxBondDim, var);

    dataFile << en1 << " " << p1 << " " << var << " " << maxBondDim << " ";

    dataFile.close();

    return 0;
}// runPT

