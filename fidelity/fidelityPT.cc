#include "itensor/all.h"

using namespace itensor;

void runPT(int, int, double, double);

int main(int argc, char *argv[])
    {
    std::vector<int> Ly={3, 5, 7, 9};
    std::vector<int> Lx={8, 16, 24, 32, 40, 48, 56, 64};
    std::vector<double> h = {2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2};
    std::vector<double> dh = {0.025, 0.01};
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
    dataFile << "en0" << " " << "M0" << " " << "p0" << " " << "var0" << " " << "maxBondDim0" << " "
             << "en1" << " " << "M1" << " " << "p1" << " " << "var1" << " " << "maxBondDim1" << " "
             << "en2" << " " << "M2" << " " << "p2" << " " << "var2" << " " << "maxBondDim2" << " "
             << "en3" << " " << "M3" << " " << "p3" << " " << "var3" << " " << "maxBondDim3" << " " << "F03" << " " << "F13" << " "
             << "en4" << " " << "M4" << " " << "p4" << " " << "var4" << " " << "maxBondDim4" << " " << "F04" << " " << "F14" << " "
             << "en5" << " " << "M5" << " " << "p5" << " " << "var5" << " " << "maxBondDim5" << " " << std::endl;


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
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, 256, 512;
    sweeps.cutoff() = 1E-10;

    //
    //solve for ground state of h-dh
    //
    // first h-dh state
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
    println("\nfirst h-dh state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en0, M0, p0, maxBondDim, var);

    dataFile << en0 << " " << M0 << " " << p0 << " " << var << " " << maxBondDim << " ";

    //
    // second h-dh state. Is the excited state for PM regime, second GS for FM regime
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
    println("\nsecond h-dh state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en1, M1, p1, maxBondDim, var);

    dataFile << en1 << " " << M1 << " " << p1 << " " << var << " " << maxBondDim << " ";

    //
    // third h-dh state. Calculate in the ferromagnetic regime to find the gap
    //
    if(true){
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
        println("\nthird h-dh state");
        printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en2, M2, p2, maxBondDim, var);

        dataFile << en2 << " " << M2 << " " << p2 << " " << var << " " << maxBondDim << " ";
        
    }

    //
    // calculate ground state of h+dh
    //
    // update Hamiltonian to correct transverse field value by adding 2*dh
    for(auto j : range1(N)){
        ampo += -4.0*dh, "Sx", j;
    }
    H = toMPO(ampo);
    auto [en3,psi3] = dmrg(H,randomMPS(sites),sweeps,{"Silent=",true});
    var = inner(psi3,H,H,psi3)-en3*en3;
    maxBondDim = maxLinkDim(psi3);
    double M3 = 0.0;
    for(auto b : range1(N)){
        psi3.position(b);
        auto m = elt( dag(prime(psi3(b),"Site")) * Sz[b-1] * psi3(b) );
        M3 += m;
    }
    M3 /= double(N);
    auto p3 = inner(psi3,P,psi3);
    auto F03 = abs(inner(psi3,psi0));
    auto F13 = abs(inner(psi3,psi1));
    println("\nfirst h+dh state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en3, M3, p3, maxBondDim, var);
    printfln("Fidelity F03 = %0.5f, Fidelity F13 = %0.5f", F03, F13);

    dataFile << en3 << " " << M3 << " " << p3 << " " << var << " " << maxBondDim << " ";
    dataFile << F03 << " " << F13 << " ";

    //
    // second d+dh state
    //
    wfs = std::vector<MPS>(1);
    wfs.at(0) = psi3;
    auto [en4,psi4] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
    var = inner(psi4,H,H,psi4)-en4*en4;
    maxBondDim = maxLinkDim(psi4);
    double M4 = 0.0;
    for(auto b : range1(N)){
        psi4.position(b);
        auto m = elt( dag(prime(psi4(b),"Site")) * Sz[b-1] * psi4(b) );
        M4 += m;
    }
    M4 /= double(N);
    auto p4 = inner(psi4,P,psi4);
    auto F04 = abs(inner(psi4,psi0));
    auto F14 = abs(inner(psi4,psi1));

    println("\nsecond h+dh state");
    printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g", en4, M4, p4, maxBondDim, var);
    printfln("Fidelity F04 = %0.5f, Fidelity F14 = %0.5f", F04, F14);

    dataFile << en4 << " " << M4 << " " << p4 << " " << var << " " << maxBondDim << " ";
    dataFile << F04 << " " << F14 << " ";

    //
    // in the ferromagnetic regime, calculate third h-dh state to find the gap
    //
    if(true){
        wfs.push_back(psi4);
        auto [en5,psi5] = dmrg(H,wfs,randomMPS(sites),sweeps,{"Silent=",true,"Weight=",20.0});
        var = inner(psi5,H,H,psi5)-en5*en5;
        maxBondDim = maxLinkDim(psi5);
        double M5 = 0.0;
        for(auto b : range1(N)){
            psi5.position(b);
            auto m = elt( dag(prime(psi5(b),"Site")) * Sz[b-1] * psi5(b) );
            M5 += m;
        }
        M5 /= double(N);
        auto p5 = inner(psi5,P,psi5);
        println("\nthird h+dh state");
        printfln("Energy = %0.3f, M = %0.3f, parity = %0.3f, maxLinkDim = %d, var = %0.3g\n", en5, M5, p5, maxBondDim, var);

        dataFile << en5 << " " << M5 << " " << p5 << " " << var << " " << maxBondDim << " ";
        
    }

    dataFile.close();

    return;
}// runPT