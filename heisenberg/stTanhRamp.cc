#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);
    
    int method = 2, N; // We assume N is even and N/2 is even.
    double v, h, quenchtau, dt;
    double tanhshift = 1.0;

    char schar1[64];
    int n1 = std::sprintf(schar1, "parameters_run%d.txt",runNumber);
    std::string s1(schar1);
    std::ifstream parameter_file(s1);
    std::string parameter;
    if ( parameter_file.is_open() ) { // always check whether the file is open
        std::getline(parameter_file, parameter); //skip header line
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        method = std::stoi(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        N = std::stoi(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        v = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        h = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        quenchtau = std::stod(parameter);
    }
    parameter_file.close();
    
    printfln("N = %d, v = %0.1f, h = %0.1f, quench tau = %0.2f", N, v, h, quenchtau);

    // We will write into a file with the time-evolved energy density at all times.
    char schar2[64];
    if(method==1){
        int n2 = std::sprintf(schar2,"N_%d_v_%0.1f_h_%0.1f_qtau_%0.2f_HeisenbergSTtanh_TEBD2.dat",N,v,h,quenchtau);
    }
    else if(method==2){
        int n2 = std::sprintf(schar2,"N_%d_v_%0.1f_h_%0.1f_qtau_%0.2f_HeisenbergSTtanh_TEBD4.dat",N,v,h,quenchtau);
    }
    else{
        printfln("Not a valid method");
        return 0;
    }
    std::string s2(schar2);
    std::ofstream enerfile;
    enerfile.open(s2); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "time" << " " << "energy" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    
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
    auto Hfinal = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,100,100,200; //gradually increase states kept
    sweeps.cutoff() = 1E-10; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  0.5*sites.op("S+",b)*sites.op("S-",b+1);
        LED[b-1] += 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        LED[b-1] += 1.0*sites.op("Sz",b)*sites.op("Sz",b+1);
    }

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, v, quenchtau, tanhshift);
    for(int b=1; b<=N; b++){
        ampo += hvals[b-1],"Sz",b;
    }
    
    // Find Initial Ground State
    auto [energy,psi] = dmrg(toMPO(ampo),initState,sweeps,{"Silent=",true});
    energy = inner(psi, Hfinal, psi);
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //calculate local energy <psi|Hf(x)|psi>
    for (int b = 1; b < N; b++){
        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        localEnergy[b-1] = elt( dag(prime(ket,"Site")) * LED[b-1] * ket);
    }
    enerfile << 0.0 << " " << energy << " " << SvN << " ";
    IndexSet bonds = linkInds(psi); //get bond dimensions
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;

    // time evolution parameters. Get time accuracy of 1E-4
    if(method == 1){ //2nd order TEBD
        dt = 0.01;
        printfln("Starting TEBD2, dt = %0.2f", dt);
    }
    else{ //4th order TEBD
        dt = 0.1;
        printfln("Starting TEBD4, dt = %0.2f", dt);
    }
    Real tval = 0.0;
    double finalTime = double(N)/2.0/v + 2.0*quenchtau*tanhshift + 2.0*tanhshift;
    int nt = int(finalTime/dt)+1;
    auto args = Args("Cutoff=",1E-10,"MaxDim=",512);
    
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, SvN, maxLinkDim(psi));

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        //update magnetic field vector
        hvals = hvector(N, tval, h, v, quenchtau, tanhshift);
        
        // TEBD time update
        std::vector<BondGate> gates;
        if(method==1){ // 2nd order TEBD
            gates = makeGates(N, hvals, dt, sites, LED);
        }
        else{ //4th order TEBD
            Real delta1 =  0.414490771794376*dt;
            Real delta2 = -0.657963087177503*dt;
            auto gatesdelta1 = makeGates(N, hvals, delta1, sites, LED);
            auto gatesdelta2 = makeGates(N, hvals, delta2, sites, LED);
            gates = gatesdelta1;
            gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
            gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
            gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
            gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
        }

        // apply Trotter gates
        gateTEvol(gates,dt,dt,psi,{args,"Verbose=",false});
        psi.orthogonalize(args); //orthogonalize to minimize bond dimensions

        // calculate energy <psi|Hf|psi>
        auto energy = innerC(psi, Hfinal, psi).real();
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi|Hf(x)|psi>
        for (int b = 1; b < N; b++){
            psi.position(b);
            auto ket = psi(b)*psi(b+1);
            localEnergy[b-1] = eltC( dag(prime(ket,"Site")) * LED[b-1] * ket ).real();
        }

        enerfile << tval << " " << energy << " " << SvN << " ";
        IndexSet bonds = linkInds(psi); //get bond dimensions
        for (int j = 0; j < N-1; j++){
            enerfile << dim(bonds[j]) << " ";
        }
        for (int j = 0; j < N-1; j++){
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, SvN, maxLinkDim(psi));
    }
    
    std::cout<< std::endl << " END PROGRAM. TIME TAKEN :";
    
    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);
    
    return 0;
}

std::vector<double> hvector(int N, double tval, double h, double v, double quenchtau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N/2; b++){
        if (b%2 == 0){
            hvals[b-1] = +h*(0.5 + 0.5*tanh(double(-b+N/2)/(v*quenchtau) - (tval-tanhshift)/quenchtau ));
        }
        else{
            hvals[b-1] = -h*(0.5 + 0.5*tanh(double(-b+N/2)/(v*quenchtau) - (tval-tanhshift)/quenchtau ));
        }
    }
        
    for (int b = N/2+1; b <= N; b++){
        if (b%2 == 0){
            hvals[b-1] = +h*(0.5 + 0.5*tanh(double(b-N/2-1)/(v*quenchtau) - (tval-tanhshift)/quenchtau ));
        }
        else{
            hvals[b-1] = -h*(0.5 + 0.5*tanh(double(b-N/2-1)/(v*quenchtau) - (tval-tanhshift)/quenchtau ));
        }
    }
    return hvals;
}

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int L, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED)
    {
    std::vector<BondGate> gates; 
    //Create the gates exp(-i*tstep/2*hterm)
    for(int i=1; i<=L; i++){
        if(i<L){
            auto hterm = LED[i-1];
            hterm += h[i-1]*op(sites,"Sz",i)*op(sites,"Id",i+1);
            auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
        else{
            auto hterm = h[i-1]*op(sites,"Id",i-1)*op(sites,"Sz",i);
            auto g = BondGate(sites,i-1,i,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
    } // for i

  //Create the gates exp(-i*tstep/2*hterm) in reverse order 
  for(int i=L; i>=1; i--){
    if(i<L){
        auto hterm = LED[i-1];
        hterm += h[i-1]*op(sites,"Sz",i)*op(sites,"Id",i+1);
        auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
    }
    else{
        auto hterm = h[i-1]*op(sites,"Id",i-1)*op(sites,"Sz",i);
        auto g = BondGate(sites,i-1,i,BondGate::tReal,dt/2.,hterm);
        gates.push_back(g);
    }
  }// for i

  return gates;
  
}// makeGates

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