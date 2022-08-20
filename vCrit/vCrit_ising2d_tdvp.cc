#include "itensor/all.h"
#include "tdvp.h"

using namespace itensor;

//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[]){
    std::clock_t tStart = std::clock();

    if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto Lx = input.getInt("Lx", 16);
    auto Ly = input.getInt("Ly", 3);
    auto h = input.getReal("h", 4.0);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxDim = input.getInt("maxDim", 128);
    
    // write results to file
    char schar[64];
    int n = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_maxDim_%d_vCrit.dat",Ly,Lx,h,maxDim);
    std::string s(schar);
    std::ofstream dataFile;
    dataFile.open(s); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "enPhi" << " " << "SvN(x,t=0)" << " " << "bondDim" << " " << std::endl;

    auto N = Ly * Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
    for(auto j : range1(N)){
        ampo += -2.0*h, "Sz", j;
    }
    auto H = toMPO(ampo);

    //initial state
    auto state = InitState(sites); 
    for(auto j : range1(N)){
        state.set(j, "Up");
    }
    auto initState = MPS(state);
    

    // 2d ising model parameters
    auto sweeps = Sweeps(15);
    sweeps.maxdim() = 10, 20, 50, maxDim;
    sweeps.cutoff() = truncE;
    sweeps.noise() = 1E-7, 1E-8, 0., 1E-7, 1E-8, 1E-9, 1E-10, 0.;

    // calculate ground state
    auto [en_psi,psi] = dmrg(H,initState,sweeps,{"Silent=",true});

    // make |phi> = Sz|psi>
    int loc = (Lx/2)*Ly+1; //centered in x and on lower row in y
    psi.position(loc);
    auto newA = 2.0*sites.op("Sz",loc) * psi(loc);
    newA.noPrime();
    auto phi = psi;
    phi.set(loc, newA);
    auto en_phi = inner(phi,H,phi);

    // make vectors to store SvN
    std::vector<Real> svn(Lx-1,0.);
    for(auto i : range1(Lx-1)){
        auto b = i*Ly;
        svn[i-1] = vonNeumannS(phi, b);
    }
    auto bonds = linkInds(phi); //get bond dimensions

    printfln("\nIteration %d, time = %0.2f; phi energy = %0.f, max link dim is %d",0,0, en_phi,maxLinkDim(phi));
    // store to file
    dataFile << 0.0 << " " << en_psi << " " << en_phi << " ";
    for(int j = 0; j<Lx-1; j++){ //save svn
        dataFile << svn[j] << " ";
    }
    for (int j=0; j<N-1; j++){
        dataFile << dim(bonds[j]) << " ";
    }
    dataFile << std::endl;

    dataFile << "tval" << " " << "enPhi" << " " << "SvN(x,t)" << " " << "bondDim" << " " << std::endl;

    // time evolution parameters
    double tval = 0.0; //time
    double ttotal = 8.0;
    Real dt = 0.1; //time step
    int Nt = int(ttotal/dt); //number of time steps
    Real delta1 =  0.414490771794376*dt; // 4th order time stepping
    Real delta2 = -0.657963087177503*dt;
    int linkCheck = int(log2( double(maxDim)) );
    int numCenter = 2; // start with two-site tdvp

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxDim;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 12;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxDim;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 12;

    printfln("\nStarting fourth order 2-site TDVP, dt = %0.2f\n", dt);

    ////////////////////////////////////////////////////////////////////////////
    ///////// time evolve //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for(int n=1; n<=Nt; n++){
        tval += dt; //update time vector
        
        // TDVP sweep
        tdvp(phi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});
        tdvp(phi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"NumCenter",numCenter});
        en_phi = tdvp(phi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});        
        
        for(auto i : range1(Lx-1)){
            auto b = i*Ly;
            svn[i-1] = vonNeumannS(phi, b);
        }
        bonds = linkInds(phi); //get bond dimensions

        //write to file
        dataFile << tval << " " << en_phi << " " ;
        for(int j = 0; j<Lx-1; j++){ //save svn
            dataFile << svn[j] << " ";
        }
        for (int j=0; j<N-1; j++){
            dataFile << dim(bonds[j]) << " ";
        }
        dataFile << std::endl;

        printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(phi));

        // check if bondDim is maxed out
        if( numCenter > 1 && dim(bonds[linkCheck-1]) >= maxDim){
            printfln("link %d has bond dimension %d", linkCheck, dim(bonds[linkCheck-1]));
            printfln("Switching to 4th order 1-site TDVP");
            numCenter = 1;
        }

    }// for n

    dataFile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

}//main

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