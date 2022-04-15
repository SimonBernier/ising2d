#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double);
//calculates local energy density
double calculateLocalEnergy(int, int, MPS, std::vector<ITensor>, std::vector<double>, SiteSet);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>, std::vector<double>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int N=16, iRange = 4; // We assume N is even and N/2 is even.
    double h=0.;

    if(argc > 1)
        N = std::stoi(argv[1]);
    if(argc > 2)
        h = std::stod(argv[2]);  
    
    printfln("N = %d, h = %0.1f", N, h);

    // We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    char schar2[128];
    int n1 = std::sprintf(schar1,"N_%d_h_%0.1f_heisR3CritEn.dat", N,h);
    int n2 = std::sprintf(schar2,"N_%d_h_%0.1f_heisR3CritSSC.dat", N,h);

    std::string s1(schar1), s2(schar2);
    std::ofstream enerfile, sscfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    sscfile.open(s2); // opens the file
    if( !sscfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    enerfile << "energy" << " " << "var" << " " << "SvN" << " " << "bondDim" << " " << "localEnergy" << " " << std::endl;
    sscfile << "szsz" << " " << "spsm" << " " << std::endl;
    
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
        if (i==1)
            g[0] = 1.0;
        else
            g[i-1] = pow( 1./double(i), 3.);
    }
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for(int i = 0; i<iRange; i++){
        for (int b = 1; b < N-i; b++){
            ampo += g[i]*0.5,"S+", b, "S-", b+i+1;
            ampo += g[i]*0.5,"S-", b, "S+", b+i+1;
            ampo += g[i]*1.0,"Sz", b, "Sz", b+i+1;
        }
    }
    if(h!=0){ //on-site field
        std::vector<double> hvals = hvector(N, h); //magnetic field vector
        for(int b=1; b<=N; b++){
            ampo += hvals[b-1],"Sz",b;
        }
    }
    auto H = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(6); //number of sweeps is 6
    sweeps.maxdim() = 10,20,50,100,100,200; //gradually increase states kept
    sweeps.cutoff() = 1E-10; //desired truncation error
    sweeps.noise() = 1E-7, 1E-8, 0;

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
    
    // Find Initial Ground State
    auto [energy,psi] = dmrg(H,initState,sweeps,{"Quiet=",true});
    auto var = inner(H, psi, H, psi) - energy*energy;
    //calculate entanglement
    auto SvN = vonNeumannS(psi, N/2);
    //get bond dimensions
    auto bonds = linkInds(psi); 
    //calculate local energy <psi|Hf(x)|psi>
    for(int b=1; b<N; b++){
        localEnergy[b-1] = calculateLocalEnergy(N, b, psi, LED, g, sites);
    }
    
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        auto [szsz,spsm] = spinspin(N/2+1,b,psi,sites);
        szszcorr[b-1] = szsz;
        spsmcorr[b-1] = spsm;
    }

    //store variables to energy file
    enerfile << energy << " " << var << " " << SvN << " ";
    for (int j=0; j<N-1; j++){
        enerfile << dim(bonds[j]) << " ";
    }
    for (int j = 0; j < N-1; j++){
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;
    //store variables to spin spin correlation file
    for (int j = 0; j < N; j++){
        sscfile << szszcorr[j] << " ";
    }
    for (int j = 0; j < N; j++){
        sscfile << spsmcorr[j] << " ";
    }
    sscfile << std::endl;
   
    enerfile.close();

    std::cout<< std::endl << " END OF PROGRAM. ";
    
    std::printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
}

// returns vector of alternating values of h
std::vector<double> hvector(int N, double h)
    {
    std::vector<double> hvals(N);
    for (int b = 1; b <= N; b++){
        if (b%2 == 0){
            hvals[b-1] = +h;
        }
        else{
            hvals[b-1] = -h;
        }
    }
    return hvals;
}

//calculate local energy density and return a vector of doubles
double calculateLocalEnergy(int N, int b, MPS psi, std::vector<ITensor> LED, std::vector<double> g, SiteSet sites){
    
    double energy;

    psi.position(b);
    auto ket = psi(b)*psi(b+1);
    energy = eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();

    for(int i = 1; i<int(size(g)); i++){

        for (int j=1; j<=i; j++){// SMART switch sites for next-nearest neighbour interaction
            for (int k=i; k>=j; k--){

                int ind = b+k+j-1;

                if (ind<N){

                    psi.position(ind);
                    auto g = BondGate(sites,ind,ind+1);
                    auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
                    AA.replaceTags("Site,1","Site,0");
                    psi.svdBond(g.i1(), AA, Fromright); //svd from the right
                    psi.position(g.i1()); //restore orthogonality center to the left

                } //if 
            } // for k
        } // for j

        for (int j=0; j<=i; j++){ //smart ordering of gates

            int ind = b+2*j;

            if (ind<N){

                psi.position(ind);
                auto ket = psi(ind)*psi(ind+1);
                energy += g[i]*eltC( dag(prime(ket,"Site")) * LED[ind-1] * ket).real();

            }
        } // for j
                
        for (int j=i; j>=1; j--){// SMART switch sites for next-nearest neighbour interaction
            for (int k=j; k<=i; k++){

                int ind = b+k+j-1;

                if (ind<N){

                    psi.position(ind);
                    auto g = BondGate(sites,ind,ind+1);
                    auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
                    AA.replaceTags("Site,1","Site,0");
                    psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
                    psi.position(g.i2()); //move orthogonality center to the right

                } // if
            } //for k
        } // for j
        psi.orthogonalize({"Cutoff=",1E-10,"MaxDim=",200}); //compress MPS to make next step faster
    }//for i

    return energy;

}//calculateLocalEnergy

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int N, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED, std::vector<double> g)
    {
    
    std::vector<BondGate> gates; 

    //Create the gates exp(-i*tstep/2*hterm)
    for(int b=1; b<=N; b++){
        if(b<N){ //nearest neighbour

            auto hterm = LED[b-1];
            hterm += h[b-1]*op(sites,"Sz",b)*op(sites,"Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
        else{

            auto hterm = h[b-1]*op(sites,"Id",b-1)*op(sites,"Sz",b);
            auto g = BondGate(sites,b-1,b,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
        for(int i = 1; i<int(size(g)); i++){ //long-range

            if (b<N-i){

                for (int j=1; j<=i; j++){// SMART switch sites for next-nearest neighbour interaction
                    for (int k=i; k>=j; k--){

                        int ind = b+k+j-1;
                        if (ind<N){
                            gates.push_back( BondGate(sites,ind,ind+1) );
                        }

                    } // for k
                } // for j

                for (int j=0; j<=i; j++){ //smart ordering of gates
                    int ind = b+2*j;
                    if (ind<N){
                        auto hterm = g[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                        auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                        gates.push_back(g);
                    }
                }
                
                for (int j=i; j>=1; j--){// SMART switch sites for next-nearest neighbour interaction
                    for (int k=j; k<=i; k++){

                        int ind = b+k+j-1;
                        if (ind<N){
                            gates.push_back( BondGate(sites,ind,ind+1) );
                        }
                        
                    } // for k
                } // for j
            } //if
        } //for i
    } // for b

    //Create the gates exp(-i*tstep/2*hterm) in reverse order 
    for(int b=N; b>=1; b--){
        for(int i = int(size(g))-1; i>=1; i--){ //long-range

            if (b<N-i){

                for (int j=1; j<=i; j++){// SMART switch sites for next-nearest neighbour interaction
                    for (int k=i; k>=j; k--){

                        int ind = b+k+j-1;
                        if (ind<N){
                            gates.push_back( BondGate(sites,ind,ind+1) );
                        }

                    } // for k
                } // for j

                for (int j=0; j<=i; j++){ //smart ordering of gates
                    int ind = b+2*j;
                    if (ind<N){
                        auto hterm = g[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                        auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                        gates.push_back(g);
                    }
                } // for j
                
                for (int j=i; j>=1; j--){// SMART switch sites for next-nearest neighbour interaction
                    for (int k=j; k<=i; k++){

                        int ind = b+k+j-1;
                        if (ind<N){
                            gates.push_back( BondGate(sites,ind,ind+1) );
                        }
                        
                    } // for k
                } // for j
            } //if
        } //for i
        if(b<N){
            auto hterm = LED[b-1];
            hterm += h[b-1]*op(sites,"Sz",b)*op(sites,"Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }
        else{
            auto hterm = h[b-1]*op(sites,"Id",b-1)*op(sites,"Sz",b);
            auto g = BondGate(sites,b-1,b,BondGate::tReal,dt/2.,hterm);
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