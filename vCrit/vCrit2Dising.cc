#include "itensor/all.h"
#include "tdvp.h"

using namespace itensor;

// calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);

int main(int argc, char *argv[])
{
  	std::clock_t tStart = std::clock();

  	if(argc < 2){ 
        printfln("Usage: %s input_file",argv[0]); 
        return 0; 
    }
    auto input = InputGroup(argv[1],"input");

    auto Ly = input.getInt("Ly", 3);
	auto Lx = input.getInt("Lx", 8);
    auto h = input.getReal("h", 1.); // 1.6700 2.54617 2.7300 2.824674
    auto truncErr = input.getReal("truncE", 1E-8);
    auto maxDim = input.getInt("maxDim", 512);
	println();

  	// We will write into a file with the time-evolved energy density at all times.
    char schar1[128];
    int n1 = std::sprintf(schar1,"Ly_%d_Lx_%d_h_%0.2f_2dTFIvCrit.dat", Ly, Lx,h);

    std::string s1(schar1);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "enPhi" << " " << "svn(x)" << " " << std::endl;

  	auto N = Ly * Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity",true});
  	auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

	// autompo hamiltonian
	auto ampo = AutoMPO(sites);
	for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
	for(int j=1; j<=N; j++){
		ampo += -2.0 * h, "Sz", j;
	}
	auto H = toMPO(ampo);

	// initial state
	auto state = InitState(sites);
	for (auto j : range1(N)){
		state.set(j, "Up");
	}

	// 2d ising model parameters
	auto sweeps = Sweeps(15);
	sweeps.maxdim() = 50, 100, 200, 200, 400, 400, maxDim;
	sweeps.cutoff() = truncErr;
	sweeps.noise() = 1E-7, 1E-8, 0.0;

	// calculate ground state
	auto [enPsi, psi] = dmrg(H, MPS(state), sweeps, {"Silent=", true});

	// make |phi> = Sz|psi>
	int loc = (Lx / 4 - 1) * Ly + 1; 
	psi.position(loc);
	auto newA = 2.0 * sites.op("Sz", loc) * psi(loc);
	newA.noPrime();
	psi.set(loc, newA);
    psi.orthogonalize({"Cutoff=",truncErr,"MaxDim=",maxDim});
    auto enPhi = innerC(psi,H,psi).real(); //energy after disturbing the ground state

	// store spin-spin correlation function
	std::vector<double> svn(Lx-1);
	// calculate the spin-spin correlation using MPS * MPO * MPS methods
	for(int j = 1; j < Lx; j++){
		auto b = j*Ly;
		svn[j-1] = vonNeumannS(psi, b);
	}

	// store to file
	dataFile << 0 << " " << enPsi << " " << maxLinkDim(psi) << " " << enPhi << " ";
	for (int j = 0; j < Lx-1; j++){
		dataFile << svn[j] << " ";
	}
	dataFile << std::endl;

	printfln("time = %0.1f; phi energy = %0.3f, max link dim is %d", 0, enPhi, maxLinkDim(psi));

	// time evolution parameters.
    double tval = 0., dt = 0.1;
    double delta1 =  0.414490771794376*dt; // for 4th order TDVP
    double delta2 = -0.657963087177503*dt;
    double finalTime = Lx/4;
    int nt=int(finalTime/dt);
    int numCenter = 2; // two-site tdvp

    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxDim;
    sweeps1.cutoff() = truncErr;
    sweeps1.niter() = 15;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxDim;
    sweeps2.cutoff() = truncErr;
    sweeps2.niter() = 15;

    println("\nStarting 4th order two-site TDVP\n");

	////////////////////////////////////////////////////////////////////////////
	///////// time evolve //////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	for (int n = 1; n <= nt; n++){
    	tval += dt; // update time

		// TDVP sweep
        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"NumCenter",numCenter});
        enPhi = tdvp(psi, H,- Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});        
        
        // calculate svn
        for(int j = 1; j < Lx; j++){
			auto b = j*Ly;
			svn[j-1] = vonNeumannS(psi, b);
		}

		// store to file
		dataFile << tval << " " << enPsi << " " << maxLinkDim(psi) << " " << enPhi << " " ;
		for(int j = 0; j<Lx-1; j++){ //save svn
			dataFile << svn[j] << " ";
		}
		dataFile << std::endl;

      	printfln("Iteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d\n", n, tval, enPhi, maxLinkDim(psi));

  	} // for n

	dataFile.close();

    print(" END OF PROGRAM. ");
    printf("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

} // main

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