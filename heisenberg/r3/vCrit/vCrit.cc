#include "itensor/all.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, double);
//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, std::vector<double>, double, SiteSet, std::vector<ITensor>, std::vector<double>);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    int N=16, maxB=512, iRange = 4; // We assume N is even and N/2 is even.
    double h=0., truncE=1E-10, alpha = 3.;

    if(argc > 1)
        N = std::stoi(argv[1]);
    if(argc > 2)
        h = std::stod(argv[2]);  
    
    printfln("N = %d, h = %0.2f", N, h);

    // We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    int n1 = std::sprintf(schar1,"N_%d_h_%0.2f_heisR3vCrit.dat", N,h);

    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "MaxDimPsi" << " " << "enPhi" << " " << "MaxDimPhi" << " " 
             << "Sz(x,y)Sz(N/2+1,0)" << " " << std::endl;
    
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
    if(h!=0){ //on-site field
        std::vector<double> hvals = hvector(N, h); //magnetic field vector
        for(int b=1; b<=N; b++){
            ampo += hvals[b-1],"Sz",b;
        }
    }
    auto H = toMPO(ampo);
    
    //sweeps
    auto sweeps = Sweeps(5); //number of sweeps is 5
    sweeps.maxdim() = 10,20,50,100,200; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error
    sweeps.noise() = 1E-7, 1E-8, 0;

    // Create the Local Energy Density Tensors
    std::vector<ITensor> LED(N-1);
    for (int b = 1; b < N; b++){
        LED[b-1] =  0.5*sites.op("S+",b)*sites.op("S-",b+1);
        LED[b-1] += 0.5*sites.op("S-",b)*sites.op("S+",b+1);
        LED[b-1] += 1.0*sites.op("Sz",b)*sites.op("Sz",b+1);
    }

    //make vector of Sz MPOs
    std::vector<MPO> Sz(N);
    for(auto j : range1(N)){
      auto ampoSz = AutoMPO(sites);
      for(auto i : range1(N)){
        if(i==j){
            ampoSz += 2.0, "Sz", i; //Sz only at site j
        }
        else{
          ampoSz += "Id", i; //identities everywhere else
        }
      }
      Sz[j-1] = toMPO(ampoSz);
    }

    // Create the SzSz and S+S- correlation vector
    std::vector<Complex> szszcorr(N,0.);
    
    // Find Initial Ground State
    auto [en_psi,psi] = dmrg(H,initState,sweeps,{"Silent=",true});

    // make |phi> = sigma_z|psi>
    int loc = N/2+1; //centered in the middle of the chain
    psi.position(loc);
    auto newA = 2.0*sites.op("Sz",loc)*psi(loc);
    newA.noPrime();
    auto phi = psi;
    phi.set(loc, newA);
    phi.orthogonalize({"Cutoff=",truncE,"MaxDim=",maxB});
    auto en_phi = innerC(phi,H,phi).real();

    //calculate the spin-spin correlation using MPS * MPO * MPS methods
    for(auto j : range1(N)){
      szszcorr[j-1] = innerC(psi, Sz[j-1], phi);
    }

    printfln("\nt=0; phi energy = %0.3f, max link dim is %d", en_phi, maxLinkDim(phi));
    // store to file
    dataFile << 0.0 << " " << en_psi << " " << maxLinkDim(psi) << " " << en_phi << " " << maxLinkDim(phi) << " ";
    for(int j = 0; j<N; j++){ //save local energycorrelation values
      dataFile << real(szszcorr[j]) << " " << imag(szszcorr[j]) << " ";
    }
    dataFile << std::endl;

    //make header for t>0 calculations
    dataFile << "t" << " " << "enPhi" << " " << "MaxDimPhi" << " " << "Sz(x,y)Sz(Lx/2,0)" << " " << std::endl;
 
    // time evolution parameters. Get time accuracy of 1E-4
    double tval = 0.0, tstep = 0.1, dt = 0.1;
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 8.;
    int nt=int(finalTime/tstep);
 
    //update magnetic field vector
    auto hvals = hvector(N, h);

    std::vector<BondGate> gates;
    auto gatesdelta1 = makeGates(N, hvals, delta1, sites, LED, J);
    auto gatesdelta2 = makeGates(N, hvals, delta2, sites, LED, J);
    gates = gatesdelta1;
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
    gates.insert(std::end(gates), std::begin(gatesdelta2), std::end(gatesdelta2));
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
    gates.insert(std::end(gates), std::begin(gatesdelta1), std::end(gatesdelta1));
    auto args = Args("Cutoff=",truncE,"MaxDim=",maxB);
    
    ////////////////////
    // TIME EVOLUTION //
    ////////////////////
    for (int n = 1; n <= nt; ++n){
        tval += tstep;

        // apply Trotter gates
        gateTEvol(gates,tstep,dt,phi,{args,"Verbose=",false});

        // calculate energy of state phi
        en_phi = innerC(phi, H, phi).real();
        
        //calculate the spin-spin correlation using MPS * MPO * MPS methods
        for(auto j : range1(N)){
            szszcorr[j-1] = innerC(exp(1_i*en_psi*tval)*psi, Sz[j-1], phi);
        }

        //write to file
        dataFile << tval << " " << en_phi << " " << maxLinkDim(phi) << " ";
        for(int j = 0; j<N; j++){ //save local energy values
            dataFile << real(szszcorr[j]) << " " << imag(szszcorr[j]) << " ";
        }
        dataFile << std::endl;

        printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(phi));

    }// for n
    
    
    dataFile.close();

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
