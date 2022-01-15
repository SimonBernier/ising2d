#include "itensor/all.h"

using namespace itensor;

void runPT(int, int, double, double);

int main(int argc, char *argv[])
    {
    std::vector<int> Ly={3, 5, 7, 9};
    std::vector<int> Lx={8, 16, 24, 32, 40, 48, 56, 64};
    std::vector<double> h = {2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2};
    std::vector<double> dh = {0.05};
    int A = Ly.size(), B = Lx.size(), C = h.size(), D = dh.size();
    int runs = A*B*C*D;
    std::vector<int> Ly_list(runs), Lx_list(runs);
    std::vector<double> h_list(runs), dh_list(runs);

    for(int a=0; a<A; a++){
        for(int b=0; b<B; b++){
            for(int c=0; c<C; c++){
                for(int d=0; d<D; d++){
                int index = a*B*C*D + b*C*D + c*D + d;
                Ly_list[index] = Ly[a];
                Lx_list[index] = Lx[b];
                h_list[index]  = h[c];
                dh_list[index] = dh[d];
                }
            }
        }
    }

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);

    runPT(Ly_list[runNumber],Lx_list[runNumber],h_list[runNumber],dh_list[runNumber]);
    return 0;
    } //main

void runPT(int Ly, int Lx, double h, double dh)
    {
    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.3f_dh_%0.4f_2dTFI_fidelity.dat",Ly,Lx,h,dh);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "en0" << " " << "M0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "M1" << " " << "var1" << " " << "maxBondDim1" << " " << "overlap01" << " "
             << "en2" << " " << "M2" << " " << "var2" << " " << "maxBondDim2" << " " << "overlap02" << " " << std::endl;


    auto N = Lx * Ly;
    auto sites = SpinHalf(N,{"ConserveSz=",false});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});
    
    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sz", j.s1, "Sz", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*(h-dh), "Sx", j;
    }

    // make vector of Sz operators
    std::vector<ITensor> Sz(N);
    for(auto j : range1(N)){
        Sz[j-1] = 2.0*sites.op("Sz",j);
    }

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, 256, 512;
    sweeps.cutoff() = 1E-12;

    //
    //solve for ground state of h-dh
    //
    auto H = toMPO(ampo);
    auto [en0,psi0] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
    PrintData(totalQN(psi0));
    auto var = inner(psi0,H,H,psi0)-en0*en0;
    auto maxBondDim = maxLinkDim(psi0);
    double M0 = 0.0;
    for(auto b : range1(N)){
        psi0.position(b);
        auto m = elt( dag(prime(psi0(b),"Site")) * Sz[b-1] * psi0(b) );
        M0 += m;
    }
    M0 /= double(N);
    auto sgnM0 = (M0 > 0) - (M0 < 0);
    println("h-dh state");
    printfln("Energy = %0.3f, FM OP = %0.3f, maxLinkDim = %d, var = %0.3g", en0, M0, maxBondDim, var);

    dataFile << en0 << " " << M0 << " " << var << " " << maxBondDim << " ";

    //
    // calculate ground state of h+dh
    //
    // update Hamiltonian to correct transverse field value by adding 2*dh
    for(auto j : range1(N)){
        ampo += -4.0*dh, "Sx", j;
    }
    H=toMPO(ampo);
    auto [en1,psi1] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
    var = inner(psi1,H,H,psi1)-en1*en1;
    maxBondDim = maxLinkDim(psi1);
    double M1 = 0.0;
    for(auto b : range1(N)){
        psi1.position(b);
        auto m = elt( dag(prime(psi1(b),"Site")) * Sz[b-1] * psi1(b) );
        M1 += m;
    }
    M1 /= double(N);
    auto sgnM1 = (M1 > 0) - (M1 < 0);
    auto F = abs(inner(psi1, psi0));
    println("first h+dh state");
    printfln("Energy = %0.3f, FM OP = %0.3f, maxLinkDim = %d, var = %0.3g", en1, M1, maxBondDim, var);
    printfln("fidelity = %0.5g", F);

    dataFile << en1 << " " << M1 << " " << var << " " << maxBondDim << " " << F << " ";
    
    //check if the algorithm finished in the same order parameter sector of the ferromagnetic regime
    //if not, find the first excited state which will return the second degenerate state
    if( (abs(M0) > 1E-4) && (sgnM0!=sgnM1) ){
        auto wfs = std::vector<MPS>(1);
        wfs.at(0) = psi1;
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
        auto sgnM2 = (M2 > 0) - (M2 < 0);
        F = abs(inner(psi2, psi0));
        println("second h+dh state");
        printfln("Energy = %0.3f, FM OP = %0.3f, maxLinkDim = %d, var = %0.3g", en2, M2, maxBondDim, var);
        printfln("overlap = %0.5g", F);

        dataFile << en2 << " " << M2 << " " << var << " " << maxBondDim << " " << F << " " << std::endl;

    }

    dataFile.close();

    return;
}// runPT