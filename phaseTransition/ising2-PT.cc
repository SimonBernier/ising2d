#include "itensor/all.h"

using namespace itensor;

void runPT(int, int);

int main(int argc, char *argv[])
    {
    std::vector<int> Ly={3, 5, 7, 9};
    std::vector<int> Lx={16, 24, 32, 48, 64};
    int A = Ly.size(), B = Lx.size();
    int runs = A*B;
    std::vector<int> Ly_list(runs), Lx_list(runs);

    for(int a=0; a<A; a++){
        for(int b=0; b<B; b++){
            int index = a*B + b;
            Ly_list[index] = Ly[a];
            Lx_list[index] = Lx[b];
        }
    }

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);

    runPT(Ly_list[runNumber],Lx_list[runNumber]);
    return 0;
    } //main

void runPT(int Ly, int Lx)
    {
    //write results to file
    char schar1[64];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_ising2dPT.dat",Ly,Lx);
    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    dataFile << "hval" << " " << "energy" << " " << "mag" << " " << "mag2" << " " << "mag4" << " " 
                << "var" << " " << "maxBondDim" << " " << std::endl;

    auto N = Lx * Ly;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // create vectors of h
    std::vector<double> h;
    if(Ly==3){
        h = {2.6, 2.62, 2.64, 2.65, 2.66, 2.665, 2.67, 2.675, 2.68, 2.685, 2.69, 2.695, 2.7, 2.8, 2.9}; //Ly=3
    }
    else if(Ly==5){
        h = {2.8, 2.84, 2.85, 2.86, 2.865, 2.87, 2.875, 2.88, 2.885, 2.89, 2.895, 2.9, 3.0, 3.1}; //Ly=5
    }
    else if(Ly==7){
        h = {2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.965, 2.97, 2.975, 2.98, 3.0, 3.1}; //Ly=7
    }
    else if(Ly==9){
        h = {2.9, 2.92, 2.94, 2.96, 2.965, 2.97, 2.975, 2.98, 2.985, 2.99, 2.995, 3.0, 3.1}; //Ly=9
    }
    else{
        println("Ly is studied only up to 9 sites");
        exit(1);
    }

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
    sweeps.maxdim() = 20, 50, 200, 400, 800;
    sweeps.cutoff() = 1E-10;

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
            mag4 += m*m*m*m;
        }
        mag /= double(N); mag2 /= double(N); mag4 /= double(N);
        printfln("Energy = %0.3f, FM OP = %0.3f, maxLinkDim = %d, var = %0.3g", en, mag, maxBondDim, var);

        dataFile << h[i] << " " << en << " " << mag << " " << mag2 << " " << mag4 << " " << var << " " << maxBondDim << " " << std::endl;

    }// for(i)

    dataFile.close();

    return;
}// runPT