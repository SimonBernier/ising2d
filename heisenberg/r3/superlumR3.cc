#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double);
//calculates local energy density
std::vector<double> calculateLocalEnergy(int, MPS, std::vector<ITensor>, SiteSet);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);
    
    int N, maxB=512, iRange = 4; // We assume N is even and N/2 is even.
    double v, h, tau, dt, truncE=1E-10;
    double tanhshift = 2.0;

    char schar1[64];
    int n1 = std::sprintf(schar1, "parameters_run%d.txt",runNumber);
    std::string s1(schar1);
    std::ifstream parameter_file(s1);
    std::string parameter;
    if ( parameter_file.is_open() ) { // always check whether the file is open
        std::getline(parameter_file, parameter); //skip header line
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        N = std::stoi(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        v = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        h = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        tau = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        truncE = std::stod(parameter);
        std::getline(parameter_file, parameter, ' '); // pipe file's content into stream
        maxB = std::stod(parameter);
    }
    parameter_file.close();
    
    printfln("N = %d, v = %0.1f, h = %0.1f, tau = %0.2f, cutoff = %0.1e, max bond dim = %d", 
                                                            N, v, h, tau, truncE, maxB);

    // We will write into a file with the time-evolved energy density at all times.
    char schar2[128];
    char schar3[128];
    int n2 = std::sprintf(schar2,"N_%d_v_%0.1f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_heisR3SuperEn.dat"
                                    ,N,v,h,tau,truncE,maxB);
    int n3 = std::sprintf(schar3,"N_%d_v_%0.1f_h_%0.1f_tau_%0.2f_cutoff_%0.0e_maxDim_%d_heisR3SuperSSC.dat"
                                    ,N,v,h,tau,truncE,maxB);

    std::string s2(schar2), s3(schar3);
    std::ofstream enerfile, sscfile;
    enerfile.open(s2); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    sscfile.open(s3); // opens the file
    if( !sscfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "time" << " " << "energy" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    sscfile << "time" << " " << "szsz" << " " << "spsm" << " " << std::endl;
    
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

    // make vector for interaction strenght depending on range
    std::vector<double> g(iRange);
    for (int i = 1; i<=iRange; i++){
        g[i-1] = pow( 1./double(i), 3.);
        printfln("%0.5f", g[i-1]);
    }
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for(int i = 1; i<=iRange; i++){
        for (int b = 1; b < N; b++){
            ampo += g[i-1]*0.5,"S+", b, "S-", b+i;
            ampo += g[i-1]*0.5,"S-", b, "S+", b+i;
            ampo += g[i-1]*1.0,"Sz", b, "Sz", b+i;
        }
    }
    auto Hfinal = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,100,200,maxB; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error

    // Create the Local Energy Density Tensors
    std::vector<double> localEnergy(N-1);
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  0.5*sites.op("S+",b)*sites.op("S-",b+1);
        LED[b-1] += 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        LED[b-1] += 1.0*sites.op("Sz",b)*sites.op("Sz",b+1);
    }

    // Create the SzSz and S+S- correlation vector
    std::vector<double> szszcorr(N), spsmcorr(N);

    //magnetic field vector
    std::vector<double> hvals = hvector(N, 0.0, h, v, tau, tanhshift);
    for(int b=1; b<=N; b++){
        ampo += hvals[b-1],"Sz",b;
    }
    
    // Find Initial Ground State
    auto [energy,psi] = dmrg(toMPO(ampo),initState,sweeps,{"Silent=",true});
    energy = inner(psi, Hfinal, psi);
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //get bond dimensions
    IndexSet bonds = linkInds(psi); 
    //calculate local energy <psi|Hf(x)|psi>
    localEnergy = calculateLocalEnergy(N, psi, LED, sites);
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        auto [szsz,spsm] = spinspin(N/2+1,b,psi,sites);
        szszcorr[b-1] = szsz;
        spsmcorr[b-1] = spsm;
    }

    //store variables to energy file
    enerfile << 0.0 << " " << energy << " " << SvN << " ";
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;
    //store variables to spin spin correlation file
    sscfile << 0.0 << " ";
    for (int j = 0; j < N; j++){
        sscfile << szszcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        sscfile << spsmcorr[j] << " ";
    }
    sscfile << std::endl;

    // time evolution parameters. Get time accuracy of 1E-4
    dt = 0.1;
    Real delta1 =  0.414490771794376*dt;
    Real delta2 = -0.657963087177503*dt;
    Real tval = 0.0;
    double finalTime = double(N)/2.0/v + 2.0*tau*(1.0+tanhshift);
    int nt = int(finalTime/dt)+1;
    auto args = Args("Cutoff=",truncE,"MaxDim=",maxB);
    
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, SvN, maxLinkDim(psi));

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += dt;

        //update magnetic field vector
        hvals = hvector(N, tval, h, v, tau, tanhshift);
        
        // 4th order TEBD time update
        std::vector<BondGate> gates;
        auto gatesdelta1 = makeGates(N, hvals, delta1, sites, LED);
        auto gatesdelta2 = makeGates(N, hvals, delta2, sites, LED);
        gates = gatesdelta1;
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
        gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
        gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));

        // apply Trotter gates
        gateTEvol(gates,dt,dt,psi,{args,"Verbose=",false});
        psi.orthogonalize(args); //orthogonalize to minimize bond dimensions
        
        // calculate energy <psi|Hf|psi>
        auto en = innerC(psi, Hfinal, psi).real();
        //calculate entanglement entropy
        SvN = vonNeumannS(psi, N/2);
        //calculate local energy <psi|Hf(x)|psi>
        localEnergy = calculateLocalEnergy(N, psi, LED, sites);

        enerfile << tval << " " << en << " " << SvN << " ";
        IndexSet bonds = linkInds(psi); //get bond dimensions
        for (int j = 0; j < N-1; j++){
            enerfile << dim(bonds[j]) << " ";
        }
        for (int j = 0; j < N-1; j++){
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, SvN, maxLinkDim(psi));

        if( n % int(1.0/dt) == 0){
            //calculate spin-spin correlation
            for (int b = 1; b <= N; b++){
                auto [szsz,spsm] = spinspin(N/2+1,b,psi,sites);
                szszcorr[b-1] = szsz;
                spsmcorr[b-1] = spsm;
            }

            //store variables to spin spin correlation file
            sscfile << tval << " ";
            for (int j = 0; j < N; j++){
                sscfile << szszcorr[j] << " ";
            }
            for (int j = 0; j < N; j++){
                sscfile << spsmcorr[j] << " ";
            }
            sscfile << std::endl;

        }//if spinspin

    }// for n
    
    enerfile.close();

    std::cout<< std::endl << " END PROGRAM. TIME TAKEN :";
    
    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

std::vector<double> hvector(int N, double tval, double h, double v, double tau, double tanhshift)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N/2; b++){
        double f = -(double(b-N/2)-0.5)/v - tval;
        if (b%2 == 0){
            hvals[b-1] = +h*(0.5 + 0.5*tanh( f/tau + tanhshift ));
        }
        else{
            hvals[b-1] = -h*(0.5 + 0.5*tanh( f/tau + tanhshift ));
        }
    }
    
    for (int b = N/2+1; b <= N; b++){
        double f = (double(b-N/2)-0.5)/v - tval;
        if (b%2 == 0){
            hvals[b-1] = +h*(0.5 + 0.5*tanh( f/tau + tanhshift ));
        }
        else{
            hvals[b-1] = -h*(0.5 + 0.5*tanh( f/tau + tanhshift ));
        }
    }
    return hvals;
}

//calculate local energy density and return a vector of doubles
std::vector<double> calculateLocalEnergy(int N, MPS psi, std::vector<ITensor> LED, SiteSet sites){
    
    const double g2 = 0.015625;
    std::vector<double> energyVector(N-1);

    for (int b = 1; b < N; b++){

        psi.position(b);
        auto ket = psi(b)*psi(b+1);
        energyVector[b-1] = eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();

        if (b<N-1){
            //switch sites for next-nearest neighbour interaction
            psi.position(b+1);
            auto g = BondGate(sites,b+1,b+2);
            auto AA = psi(b+1)*psi(b+2)*g.gate(); //contract over bond b+1
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
            psi.position(g.i1()); //restore orthogonality center to the left 

            ket = psi(b)*psi(b+1);
            energyVector[b-1] += g2*eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();

            //switch back sites
            AA = psi(b+1)*psi(b+2)*g.gate(); //contract over bond b+1
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
            psi.position(g.i2()); //move orthogonality center to the right
        }

    }//for b

    return energyVector;

}//calculateLocalEnergy

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int L, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED)
    {
    
    const double g2 = 0.015625;
    std::vector<BondGate> gates; 
    //Create the gates exp(-i*tstep/2*hterm)
    for(int i=1; i<=L; i++){
        if(i<L){ //nearest neighbour
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
        if(i<L-1 && i%2!=0){ //next-nearest neighbour
            gates.push_back( BondGate(sites, i+1, i+2) ); //swap sites
            
            auto hterm = LED[i-1]; //time evolve sites i, i+2
            auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,g2*hterm);
            gates.push_back(g);
            
            hterm = LED[i+1]; //time evolve site i+1, i+3
            g = BondGate(sites,i+2,i+3,BondGate::tReal,dt/2.,g2*hterm);
            gates.push_back(g);

            gates.push_back( BondGate(sites, i+1, i+2) ); //swap back
        }
    } // for i

    //Create the gates exp(-i*tstep/2*hterm) in reverse order 
    for(int i=L; i>=1; i--){
        if(i<L-1 && i%2!=0){ //next-nearest neighbour
            gates.push_back( BondGate(sites, i+1, i+2) ); //swap sites
            
            auto hterm = LED[i-1]; //time evolve
            auto g = BondGate(sites,i,i+1,BondGate::tReal,dt/2.,g2*hterm);
            gates.push_back(g);

            hterm = LED[i+1]; //time evolve site i+1, i+3
            g = BondGate(sites,i+2,i+3,BondGate::tReal,dt/2.,g2*hterm);
            gates.push_back(g);

            gates.push_back( BondGate(sites, i+1, i+2) ); //swap back
        }
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
    }// for i, nearest-neighbour 

  return gates;
  
}// makeGates

//calculate entanglement
Real vonNeumannS(MPS psi, int b){
    Real SvN = 0.;

    //choose orthogonality center and perform svd
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

//calculate spin-spin correlator
std::tuple<double, double> spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrZ, corrPM;

    psi.position(b);
    if(b>center){ //bring site b next to the center from right
        for(int n=b-1; n>center; n--){
            auto g = BondGate(sites,n,n+1);
            auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
            AA.replaceTags("Site,1","Site,0");
            psi.svdBond(g.i1(), AA, Fromright); //svd from the right
            psi.position(g.i1()); //move orthogonality center to the left 
        }
        auto ket = psi(center)*psi(center+1);
        auto SzSz = sites.op("Sz",center)*sites.op("Sz",center+1);
        auto SpSm = sites.op("S+",center+1)*sites.op("S-",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }
    else if(b<center){ //bring site b next to the center from left
        for(int n=b; n<center-1; n++){
          auto g = BondGate(sites,n,n+1);
          auto AA = psi(n)*psi(n+1)*g.gate(); //contract over bond n
          AA.replaceTags("Site,1","Site,0");
          psi.svdBond(g.i1(), AA, Fromleft); //svd from the right
          psi.position(g.i2()); //move orthogonality center to the right 
        }
        auto ket = psi(center-1)*psi(center);
        auto SzSz = sites.op("Sz",center-1)*sites.op("Sz",center);
        auto SpSm = sites.op("S+",center-1)*sites.op("S-",center);
        corrZ = eltC( dag(prime(ket,"Site")) * SzSz * ket).real();
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }
    else{
        corrZ = 0.25;
        auto ket = psi(center);
        auto SpSm = 0.5*sites.op("Id",center) + sites.op("Sz",center);
        corrPM = eltC( dag(prime(ket,"Site")) * SpSm * ket).real();
    }

    return {corrZ, corrPM};

}//Szsz