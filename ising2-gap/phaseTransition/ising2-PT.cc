#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[])
    {
    int Nx = 64;
    int Ny = 3;
    if(argc > 2)
        Ny = std::stoi(argv[2]);
    if(argc > 1)
        Nx = std::stoi(argv[1]);

    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Nx_%d_Ny_%d_ising2dPT.dat",Nx,Ny);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    auto N = Nx * Ny;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Nx, Ny, {"YPeriodic = ", true});

    // create vectors of h
    std::vector<double> h = {2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2};
    int iter = h.size();
    std::vector<double> diff(iter,0.0);
    for(int i=1; i<iter; i++){
        diff[i] = h[i]-h[i-1];
    }
    
    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sz", j.s1, "Sz", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*h[0], "Sx", j;
    }

    // make vector of Sz operators
    std::vector<ITensor> Sz(N);
    for(auto j : range1(N)){
        Sz[j-1] = 2.0*sites.op("Sz",j);
    }

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, 200, 400;
    sweeps.cutoff() = 1E-10;

    dataFile << "hval" << " " << "energy" << " " << "mag" << " " << "mag2" << " " << "mag4" << " " << "var" << " " << "maxBondDim" << " " << std::endl;
    
    //
    // loop over values of h
    //
    for(int i=0; i<iter; i++){

        printfln("\nStarting Iteration %d",i+1);

        for(auto j : range1(N)){
            ampo += -2.0*diff[i], "Sx", j;
        }

        auto H = toMPO(ampo); //12x12 matrices
        auto [en,psi] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
        auto var = inner(psi,H,H,psi)-en*en;
        auto maxBondDim = maxLinkDim(psi);
        double mag = 0.0, mag2 = 0.0, mag4 = 0.0;
        for(auto b : range1(N)){
            psi.position(b);
            auto m = elt( dag(prime(psi(b),"Site")) * Sz[b-1] * psi(b) );
            mag += m;
            mag2 += m*m;
            mag4 = m*m*m*m;
        }
        mag /= double(N); mag2 /= double(N); mag4 /= double(N);
        printfln("Energy = %0.3f, magnetization = %0.3f, maxLinkDim = %d, var = %0.3g", en, mag, maxBondDim, var);

        dataFile << h[i] << " " << en << " " << mag << " " << mag2 << " " << mag4 << " " << var << " " << maxBondDim << " " << std::endl;

    }// for(i)

    dataFile.close();

    return 0;
}// main
