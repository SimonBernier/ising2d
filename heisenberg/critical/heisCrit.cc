#include "itensor/all.h"

using namespace itensor;

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int N = 16;
    if(argc > 1)
        N = std::stoi(argv[1]);
       
    printfln("N = %d", N);

    // We will write into a file with the time-evolved energy density at all times.
    char schar[64];
    int n1 = std::sprintf(schar,"N_%d_HeisenbergCritical.dat",N);
    std::string s1(schar);
    std::ofstream enerfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "energy" << " " << "var" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    
    auto sites = SpinHalf(N);
    auto state = InitState(sites);
    for(int i = 1; i <= N; i++){
        if(i%2 == 0)
            state.set(i,"Dn");
        else
            state.set(i,"Up");
    }
    auto initState = MPS(state);
    PrintData(totalQN(initState));
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for (int b = 1; b < N; b++){
        ampo += 0.5,"S+", b, "S-", b+1;
        ampo += 0.5,"S-", b, "S+", b+1;
        ampo += 1.0,"Sz", b, "Sz", b+1;
    }
    auto H = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(10); //number of sweeps is 5
    sweeps.maxdim() = 10,20,100,100,200,200,512; //gradually increase states kept
    sweeps.cutoff() = 1E-10; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  0.5*sites.op("S+",b)*sites.op("S-",b+1);
        LED[b-1] += 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        LED[b-1] += 1.0*sites.op("Sz",b)*sites.op("Sz",b+1);
    }
    
    // Find Initial Ground State
    auto [energy,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto var = inner(H, psi, H, psi) - energy*energy;
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|H(x)|psi>
    for (int b = 1; b < N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        localEnergy[b-1] = elt( dag(prime(ket,"Site")) * LED[b-1] * ket);
    }
    enerfile << energy << " " << var << " " << SvN << " ";
    IndexSet bonds = linkInds(psi); //get bond dimensions
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;

    printfln("N = %d, energy = %0.3f, SvN = %0.3f, maxDim = %d", N, energy, SvN, maxLinkDim(psi));

    std::cout<< std::endl << " END PROGRAM. TIME TAKEN :";
    
    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);
    
    return 0;
}

Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //calculate entanglement
    psi.position(b);
    auto l = leftLinkIndex(psi,b);
    auto s = siteIndex(psi,b);
    auto [U,S,V] = svd(psi(b),{l,s});
    auto u = commonIndex(U,S);

    //Apply von Neumann formula
    //to the squares of the singular values
    for(auto n : range1(dim(u))){
        auto Sn = elt(S,n,n);
        auto p = sqr(Sn);
        if(p > 1E-12) SvN += -p*log(p);
    }
    return SvN;

}//vonNeumannS