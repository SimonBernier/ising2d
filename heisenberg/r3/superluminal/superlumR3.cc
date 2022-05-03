#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double, double, double, double, double);
//calculates local energy density
double calculateLocalEnergy(int, int, MPS, std::vector<ITensor>, std::vector<double>, SiteSet);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>, std::vector<double>);
// update time dependent gates
void updateGates(int, std::vector<BondGate>&, std::vector<double>, double, SiteSet, std::vector<ITensor>, std::vector<double>);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
std::tuple<double, double, double> spinspin(int,int,MPS,SiteSet);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int runNumber = 0;
    if(argc > 1)
        runNumber = std::stoi(argv[1]);
    
    int N, maxB=512, iRange = 4; // We assume N is even and N/2 is even.
    double v, h, tau, truncE=1E-10, alpha = 3.0, tanhshift = 2.0;

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
    sscfile << "time" << " " << "szsz" << " " << "spsm" << " " << "sz" << " " << std::endl;
    
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
    std::vector<double> J(iRange);
    for (int i = 1; i<=iRange; i++){
        if (i==1)
            J[0] = 1.0;
        else
            J[i-1] = pow( 1./double(i), alpha);
    }
    
    // Create the Target Hamiltonian and find the Ground State Energy Density
    auto ampo = AutoMPO(sites);
    
    for(int i = 0; i<iRange; i++){
        for (int b = 1; b < N-i; b++){
            ampo += J[i]*0.5,"S+", b, "S-", b+i+1;
            ampo += J[i]*0.5,"S-", b, "S+", b+i+1;
            ampo += J[i]*1.0,"Sz", b, "Sz", b+i+1;
        }
    }
    auto Hfinal = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,50,100; //gradually increase states kept
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
    std::vector<double> szszcorr(N), spsmcorr(N), expSz(N);

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
    for(int b=1; b<N; b++){
        localEnergy[b-1] = calculateLocalEnergy(N, b, psi, LED, J, sites);
    }
    //calculate spin-spin correlation
    for (int b = 1; b <= N; b++){
        auto [szsz,spsm,sz] = spinspin(N/2+1,b,psi,sites);
        szszcorr[b-1] = szsz;
        spsmcorr[b-1] = spsm;
        expSz[b-1] = sz;
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
    for (int j = 0; j < N; j++){
        sscfile << expSz[j] << " ";
    }
    sscfile << std::endl;

    // time evolution parameters. Get time accuracy of 1E-4
    double tval = 0.0;
    double dt = 0.125;
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 0.1*double(N)/1.5707963 + 0.5*double(N)/v + 2.0*tau*tanhshift; // 0.1*N/c + 0.5*N/v + 2*tau*shift
    int nt = int(finalTime/dt);

    auto args = Args("Cutoff=",truncE,"MaxDim=",maxB);
    
    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, SvN, maxLinkDim(psi));

    auto gatesdelta1 = makeGates(N, hvals, delta1, sites, LED, J);
    auto gatesdelta2 = makeGates(N, hvals, delta2, sites, LED, J);

    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; n++){

        tval += dt;

        //update h field
        hvals = hvector(N, tval, h, v, tau, tanhshift);

        //make gates
        std::vector<BondGate> gates;    
        updateGates(N, gatesdelta1, hvals, delta1, sites, LED, J);
        updateGates(N, gatesdelta2, hvals, delta2, sites, LED, J);
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
        for(int b=1; b<N; b++){
            localEnergy[b-1] = calculateLocalEnergy(N, b, psi, LED, J, sites);
        }

        enerfile << tval << " " << en << " " << SvN << " ";
        IndexSet bonds = linkInds(psi); //get bond dimensions
        for (int j = 0; j < N-1; j++){
            enerfile << dim(bonds[j]) << " ";
        }
        for (int j = 0; j < N-1; j++){
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        printfln("t = %0.3f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, SvN, maxLinkDim(psi));

        if( n % int(1.0/dt) == 0){
            //calculate spin-spin correlation
            for (int b = 1; b <= N; b++){
                auto [szsz,spsm,sz] = spinspin(N/2+1,b,psi,sites);
                szszcorr[b-1] = szsz;
                spsmcorr[b-1] = spsm;
                expSz[b-1] = sz;
            }

            //store variables to spin spin correlation file
            sscfile << tval << " ";
            for (int j = 0; j < N; j++){
                sscfile << szszcorr[j] << " ";
            }
            for (int j = 0; j < N; j++){
                sscfile << spsmcorr[j] << " ";
            }
            for (int j = 0; j < N; j++){
                sscfile << expSz[j] << " ";
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
double calculateLocalEnergy(int N, int b, MPS psi, std::vector<ITensor> LED, std::vector<double> J, SiteSet sites){
    
    double energy;

    psi.position(b);
    auto ket = psi(b)*psi(b+1);
    energy = eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();

    for(int i = 1; i<int(size(J)); i++){

        if (b+i+1 <= N){

            for (int j=i; j>=1; j--){// switch sites for long-range interaction

                int ind = b+j;
                psi.position(ind);
                auto g = BondGate(sites,ind,ind+1);
                auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
                AA.replaceTags("Site,1","Site,0");
                psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
                psi.position(g.i1()); //restore orthogonality center to the left

            } // for j

            psi.position(b);
            ket = psi(b)*psi(b+1);
            energy += J[i]*eltC( dag(prime(ket,"Site")) * LED[b-1] * ket).real();
                    
            for (int j=1; j<=i; j++){// SMART switch sites for next-nearest neighbour interaction

                int ind = b+j;

                psi.position(ind);
                auto g = BondGate(sites,ind,ind+1);
                auto AA = psi(ind)*psi(ind+1)*g.gate(); //contract over bond ind
                AA.replaceTags("Site,1","Site,0");
                psi.svdBond(g.i1(), AA, Fromleft); //svd from the left
                psi.position(g.i2()); //move orthogonality center to the right

            } // for j
        }// if
    }//for i

    return energy;

}//calculateLocalEnergy

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int N, std::vector<double> h, double dt, SiteSet sites, std::vector<ITensor> LED, std::vector<double> J)
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
        }// if
    }// for b

    for(int i = 1; i<int(size(J)); i++){ //long-range
        //printfln("\ninteraction range %d", i+1);

	for( int b=1; b<N-i; b+=i+1){

	    int skip=0, nsg=i;
	    if ( (N-b+1) < 2*(i+1)){
		skip = (2*(i+1)-(N-b+1))%(i+1);
		nsg = i+1 - skip;
		//printf("(skip %d of %d)", skip, i+1);
	    }

            for (int j=1; j<=nsg; j++){// SMART switch sites for next-nearest neighbour interaction
                for (int k=i; k>=j; k--){

                    int ind = b+k+j-1;
		            //printf(" sg%d ", ind);
                    gates.push_back( BondGate(sites,ind,ind+1) );

                } // for k
            } // for j

            for (int j=0; j<=i-skip; j++){ //smart ordering of gates
                //printf(" %d-%d ", b+j, b+j+i+1);
		        int ind = b+2*j;
                auto hterm = J[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                gates.push_back(g);
            }// for j

            for (int j=nsg; j>=1; j--){// SMART switch sites for next-nearest neighbour interaction
                for (int k=j; k<=i; k++){

                    int ind = b+k+j-1;
		            //printf(" sg%d ", ind);
                    gates.push_back( BondGate(sites,ind,ind+1) );
                } // for k
            } // for j
        } //for b
    } // for i

    //Create the gates exp(-i*tstep/2*hterm) in reverse order
    for(int i = int(size(J))-1; i>=1; i--){ //long-range
        //printfln("\ninteraction range %d", i+1);

        for(int b=1+(i+1)*((N-1)/(i+1)-1); b>=1; b-=i+1){

	    int skip=0, nsg=i;
            if ( (N-b+1) < 2*(i+1)){
                skip = (2*(i+1)-(N-b+1))%(i+1);
                nsg = i+1 - skip;
		        //printf("(skip %d of %d)", skip, i+1);
            }

            for (int j=1; j<=nsg; j++){// SMART switch sites for next-nearest neighbour interaction
                for (int k=i; k>=j; k--){

                    int ind = b+k+j-1;
		            //printf(" sg%d ", ind);
                    gates.push_back( BondGate(sites,ind,ind+1) );

                } // for k
            } // for j

            for (int j=i-skip; j>=0; j--){ //smart ordering of gates

                int ind = b + 2*j;
                //printf(" %d-%d ", b+j, b+j+i+1);

		        auto hterm = J[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                gates.push_back(g);

            } // for j

            for (int j=nsg; j>=1; j--){// SMART switch sites for next-nearest neighbour interaction
                for (int k=j; k<=i; k++){

                    int ind = b+k+j-1;
                    //printf(" sg%d ", ind);
                    gates.push_back( BondGate(sites,ind,ind+1) );

                } // for k
            } // for j
        } //for b
    }// for i

    for(int b=N; b>=1; b--){ // nearest-neighbour
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
    }// for b

  return gates;

}// makeGates

// update time dependent gates
void updateGates(int N, std::vector<BondGate>& gates, std::vector<double> h, 
                                  double dt, SiteSet sites, std::vector<ITensor> LED, std::vector<double> J)
    {

    int indG=0;

    //Create the gates exp(-i*tstep/2*hterm)
    for(int b=1; b<=N; b++){
        if(b<N){ //nearest neighbour

            auto hterm = LED[b-1];
            hterm += h[b-1]*op(sites,"Sz",b)*op(sites,"Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,dt/2.,hterm);
            gates[indG] = g;
            indG++;
        }
        else{

            auto hterm = h[b-1]*op(sites,"Id",b-1)*op(sites,"Sz",b);
            auto g = BondGate(sites,b-1,b,BondGate::tReal,dt/2.,hterm);
            gates[indG] = g;
            indG++;
        }// if
    }// for b

    for(int i = 1; i<int(size(J)); i++){ //long-range

        for( int b=1; b<N-i; b+=i+1){

            int skip=0;
            if ( (N-b+1) < 2*(i+1)){
                skip = (2*(i+1)-(N-b+1))%(i+1);
            }
            indG += (i+1-skip)*(i+skip)/2;

            for (int j=0; j<=i-skip; j++){ //smart ordering of gates
		        int ind = b+2*j;
                auto hterm = J[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                gates[indG] = g;
                indG++;
            }// for j

            indG += (i+1-skip)*(i+skip)/2;

        } //for b
    } // for i

    //Create the gates exp(-i*tstep/2*hterm) in reverse order
    for(int i = int(size(J))-1; i>=1; i--){ //long-range
        //printfln("\ninteraction range %d", i+1);

        for(int b=1+(i+1)*((N-1)/(i+1)-1); b>=1; b-=i+1){

	        int skip=0;
            if ( (N-b+1) < 2*(i+1)){
                skip = (2*(i+1)-(N-b+1))%(i+1);
            }
            indG += (i+1-skip)*(i+skip)/2;

            for (int j=i-skip; j>=0; j--){ //smart ordering of gates

                int ind = b + 2*j;
		        auto hterm = J[i]*LED[ind-1]; //time evolve sites b+j, b+j+i+1
                auto g = BondGate(sites,ind,ind+1,BondGate::tReal,dt/2.,hterm);
                gates[indG] = g;
                indG++;

            } // for j

            indG += (i+1-skip)*(i+skip)/2;

        } //for b
    }// for i

    for(int b=N; b>=1; b--){ // nearest-neighbour
        if(b<N){
            auto hterm = LED[b-1];
            hterm += h[b-1]*op(sites,"Sz",b)*op(sites,"Id",b+1);
            auto g = BondGate(sites,b,b+1,BondGate::tReal,dt/2.,hterm);
            gates[indG] = g;
            indG++;
        }
        else{
            auto hterm = h[b-1]*op(sites,"Id",b-1)*op(sites,"Sz",b);
            auto g = BondGate(sites,b-1,b,BondGate::tReal,dt/2.,hterm);
            gates[indG] = g;
            indG++;
        }
    }// for b

    return;

}// updateGates

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
std::tuple<double, double, double> spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrZ, corrPM, expZ;

    psi.position(b);
    expZ = eltC(dag(prime(psi(b),"Site")) * sites.op("Sz",b) * psi(b)).real();
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

    return {corrZ, corrPM, expZ};

}//Szsz
