#include "itensor/all.h"

using namespace itensor;

void runPT(int, int);

int main(int argc, char *argv[])
    {
    std::vector<int> Ly={2, 3, 4, 5, 6};
    std::vector<int> AR={2, 3, 4, 5};
    int A = Ly.size(), B = AR.size();
    int runs = A*B;
    std::vector<int> Ly_list(runs), AR_list(runs);

    for(int a=0; a<A; a++){
        for(int b=0; b<B; b++){
            int index = a*B + b;
            Ly_list[index] = Ly[a];
            AR_list[index] = AR[b];
        }
    }

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);

    runPT(Ly_list[runNumber],AR_list[runNumber]);
    return 0;
    } //main

void runPT(int Ly, int AR)
    {
    int Lx = AR*Ly;
    double PI = 3.1415926535;
    
    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_ising2dGap.dat",Ly,Lx);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "h" << " " << "en0" << " " << "SR0" << " " << "var0" << " " << "maxBondDim0" << " "
                           << "en1" << " " << "SR1" << " " << "var1" << " " << "maxBondDim1" << " " << std::endl;

    auto N = Lx * Ly;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // create vectors of h
    std::vector<double> h = {3.0, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10};

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

    // make vector for spin-reveral operator S
    std::vector<ITensor> S(N);
    for(auto j : range1(N)){
        S[j-1] = sites.op("Sx",j);
    }

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, 200;
    sweeps.cutoff() = 1E-10;

    //
    // loop over values of h
    //
    for(int i=0; i<iter; i++){

        printfln("\nStarting Iteration %d",i+1);

        for(auto j : range1(N)){
            ampo += -2.0*diff[i], "Sx", j;
        }

        auto H = toMPO(ampo);
        auto [en0,psi0] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
        auto var0 = inner(psi0,H,H,psi0)-en0*en0;
        auto maxBondDim0 = maxLinkDim(psi0);

        //calculate spin-reversal
        double SR0 = 0.0;
        for(auto b : range1(N)){
            psi0.position(b);
            auto m = elt( dag(prime(psi0(b),"Site")) * S[b-1] * psi0(b) );
            SR0 += m;
        }
        SR0 = real(exp(0.5_i*PI*(SR0+Lx*Ly)));
        printfln("GS energy = %0.3f, SR = %0.3f, maxLinkDim = %d, var = %0.3g", en0, SR0, maxBondDim0, var0);

        //
        auto wfs = std::vector<MPS>(1);
        wfs.at(0) = psi0;

        //excited state with energy penalty to psi0
        auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
        auto var1 = inner(psi1,H,H,psi1)-en1*en1;
        auto maxBondDim1 = maxLinkDim(psi1);

        //calculate spin-reversal
        double SR1 = 0.0;
        for(auto b : range1(N)){
            psi1.position(b);
            auto m = elt( dag(prime(psi1(b),"Site")) * S[b-1] * psi1(b) );
            SR1 += m;
        }
        SR1 = real(exp(0.5_i*PI*(SR1+Lx*Ly)));
        printfln("ES energy = %0.3f, SR = %0.3f, maxLinkDim = %d, var = %0.3g", en1, SR1, maxBondDim1, var1);

        dataFile << h[i] << " " << en0 << " " << SR0 << " " << var0 << " " << maxBondDim0 << " "
                                << en1 << " " << SR1 << " " << var1 << " " << maxBondDim1 << " " << std::endl;

    }// for(i)

    dataFile.close();

    return;
}// runPT