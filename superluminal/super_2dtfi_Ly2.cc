#include "itensor/all.h"
#include "tdvp.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, int, double, double, double, double, double);
//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS, std::vector<ITensor>,
                                                                 std::vector<std::vector<ITensor>>);
// calculate energy of horizontal bonds for one column
std::vector<double> calculateLRenergy(int, int, MPS, std::vector<std::vector<ITensor>>, SiteSet);
//calculate Von Neumann entanglement entropy
Real vonNeumannS(MPS, int);
//calculate spin-spin correlator
double spinspin(int,int,MPS,SiteSet);

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
    auto hc = input.getReal("hc",2.5);
    auto tau = input.getReal("tau",1.0);
    auto v = input.getReal("v",3.0);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxDim = input.getInt("maxDim", 128);
    auto tanhshift = input.getReal("tanhshift",2.0);
    auto dt = input.getReal("dt",0.1);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_v_%0.2f_tau_%0.1f_maxDim_%d.dat",Ly,Lx,h,v,tau,maxDim);

    std::string s1(schar);
    std::ofstream dataFile;
    dataFile.open(s1); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    dataFile << "tval" << " " << "energy" << " " << "svn" << " " << "MaxDim" << " " 
                << "localEnergy" << " " << "sxsx" << " " << std::endl;

    auto N = Ly * Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
    //final Hamiltonian
    for(auto j : range1(N)){
        ampo += -2.0*hc, "Sz", j;
    }
    auto Hfinal = toMPO(ampo);

    //initial state
    auto state = InitState(sites); 
    for(auto j : range1(N)){
        state.set(j, "Up");
    }
    auto initState = MPS(state);

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, maxDim;
    sweeps.cutoff() = 1E-10;

    // calculate initial local energy density
    std::vector<double> localEnergy(Ly*(Lx-1),0.0); // local energy density vector
    std::vector<double> sxsx(N,0.0);

    //make 2D vector of ITensor for local energy operators
    //long-range interactions have the same structure as nearest-neighbour when we use swap gates
    std::vector<ITensor> LED(Lx);
    std::vector<std::vector<ITensor>> LED_LR(Lx-1, std::vector<ITensor>(Ly));
    // make local energy tensors
    for(int i=1; i<=Lx; i++){

        int b = (i-1)*Ly+1;

        // MPS nearest-neighbour
        LED[i-1] = -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);

        //MPS long-range
        if(i<Lx){
            LED_LR[i-1][0] = -4.0*sites.op("Sx",b)*sites.op("Sx",b+1);
            LED_LR[i-1][1] = -4.0*sites.op("Sx",b+Ly)*sites.op("Sx",b+Ly+1);
            LED_LR[i-1][0] += -hc*sites.op("Sz",b)*sites.op("Id",b+1);
            LED_LR[i-1][0] += -hc*sites.op("Id",b)*sites.op("Sz",b+1);
            LED_LR[i-1][1] += -hc*sites.op("Sz",b+Ly)*sites.op("Id",b+Ly+1);
            LED_LR[i-1][1] += -hc*sites.op("Id",b+Ly)*sites.op("Sz",b+Ly+1);
            if( i==1 ){ // add to left
                LED_LR[i-1][0] += -hc*sites.op("Sz",b)*sites.op("Id",b+1);
                LED_LR[i-1][1] += -hc*sites.op("Sz",b+Ly)*sites.op("Id",b+Ly+1);
            }
            else if( i==Lx-1){ // add to right
                LED_LR[i-1][0] += -hc*sites.op("Id",b)*sites.op("Sz",b+1);
                LED_LR[i-1][1] += -hc*sites.op("Id",b+Ly)*sites.op("Sz",b+Ly+1);
            }
        }
    }

    // revert to H at t=0 (gapped ground state)
    auto hvals = hvector(Lx, Ly, 0.0, h, v, tau, tanhshift);
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<=Ly; j++){
            int b = Ly*(i-1) + j;
            ampo += -2.0*hvals[b-1], "Sz", b;
        }
    }
    auto H = toMPO(ampo);

    //DMRG to find ground state at t=0
    auto [energy,psi] = dmrg(H,initState,sweeps,{"Silent=",true});
    energy = inner(psi,Hfinal,psi); //energy <psi(0)|H_final|psi(0)>
    // calculate von Neumann S
    auto svN = vonNeumannS(psi, N/2);
    // calculate local energy density
    localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LED_LR);
    for(int b = 1; b<=N; b++){
        sxsx[b-1] = spinspin(N/2+1, b, psi, sites);
    }

    //store to file
    dataFile << 0.0 << " " << energy << " " << svN << " " << maxLinkDim(psi) << " "; //print to file
    for(int j = 0; j<Ly*(Lx-1); j++){ //save local energy values
        dataFile << localEnergy[j] << " ";
    }
    for(int j = 0; j<N; j++){ //save local energy values
        dataFile << sxsx[j] << " ";
    }
    dataFile << std::endl;

    // time evolution parameters
    double tval = 0.0; //time
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 0.25*double(Lx) + 0.5*double(Lx)/v + 2.0*tau*tanhshift; // 0.5*N/2 + 0.5*N/v + 2*tau*shift
    int nt = int(finalTime/dt);
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

    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, svN, maxLinkDim(psi));

    ////////////////////////////////////////////////////////////////////////////
    ///////// time evolve //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for(int n=1; n<=nt; n++){

        tval += dt; //update time vector

        auto htemp = hvector(Lx, Ly, tval, h, v, tau, tanhshift);
        for(int i=1; i<=Lx; i++){
            for(int j=1; j<=Ly; j++){
                int b = Ly*(i-1) + j;
                auto hDiff = htemp[b-1] - hvals[b-1]; // local difference in h
                ampo += -2.0*hDiff, "Sz", b; // update MPO
                hvals[b-1] = htemp[b-1]; // change to new field for next step
            }
        }
        H = toMPO(ampo);

        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"NumCenter",numCenter});
        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"NumCenter",numCenter});        
        
        // calculate energy <psi(t)|H_final|psi(t)>
        energy = innerC(psi,Hfinal,psi).real();
        //calculate entanglement entropy
        svN = vonNeumannS(psi, N/2);
        // calculate local energy density <psi(t)|H_final(x,y)|psi(t)>
        localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LED_LR);
        for(int b = 1; b<=N; b++){
            sxsx[b-1] = spinspin(N/2+1, b, psi, sites);
        }

        //write to file
        dataFile << tval << " " << energy << " " << svN << " " << maxLinkDim(psi) << " ";
        for(int j = 0; j<Ly*(Lx-1); j++){ //save local energy values
            dataFile << localEnergy[j] << " ";
        }
        for(int j = 0; j<N; j++){ //save local energy values
            dataFile << sxsx[j] << " ";
        }
        dataFile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, svN, maxLinkDim(psi));

        // check if bondDim is maxed out
        auto bonds = linkInds(psi); //get bond dimensions
        if( numCenter > 1 && dim(bonds[linkCheck-1]) >= maxDim){
            printfln("link %d has bond dimension %d", linkCheck, dim(bonds[linkCheck-1]));
            printfln("Switching to 4th order 1-site TDVP");
            numCenter = 1;
        }
    }

    dataFile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;
    }

// transverse field vector
std::vector<double> hvector(int Lx, int Ly, double tval, double h, double v, double tau, double tanhshift){

    std::vector<double> hvals(Lx*Ly);

    for (int i = 1; i <= Lx; i++){

        double f = abs(double(i-Lx/2)-0.5)/v - tval;

        for (int j = 1; j <= Ly; j++){

            int index = Ly*(i-1) + j;
            hvals[index-1] = h*(0.5 + 0.5*tanh( f/tau + tanhshift ));

        } // for 
        
    }// for i

    return hvals;
}

// calculates local energy for 2D MPS using gates
std::vector<double> calculateLocalEnergy(int Lx, int Ly, SiteSet sites, MPS psi,
                                std::vector<ITensor> LED, 
                                std::vector<std::vector<ITensor>> LED_LR){
  
    std::vector<double> vertBonds(Lx, 0.0);
    std::vector<std::vector<double>> horBonds((Lx-1), std::vector<double>(Ly, 0.0)); 
    std::vector<double> localEnergy( (Lx-1)*Ly, 0.0); // interpolated energy density

    // MPS nearest-neighbour interaction
    for(int i=1; i<=Lx; i++){
        int index = (i-1)*Ly + 1;

        psi.position(index);
        auto ket = psi(index)*psi(index+1);
        vertBonds[i-1] = eltC(dag(prime(ket,"Site")) * LED[i-1] * ket).real();

    }// for i

    // MPS long-range interactions
    for(int i=1; i<Lx; i++){

        auto lrEnergy = calculateLRenergy(i, Ly, psi, LED_LR, sites);

        for(int m = 0; m<Ly; m++){

            horBonds[i-1][0] = lrEnergy[0];
            horBonds[i-1][1] = lrEnergy[1];

        } // for m

    }// for i
    
    // interpolate into localEnergy
    for(int i=1; i<Lx; i++){
        for(int j=1; j<=Ly; j++){
            int b = (i-1)*Ly + j;
            localEnergy[b-1] = horBonds[i-1][j-1];
            localEnergy[b-1] += 0.25 * ( vertBonds[i-1] + vertBonds[i] );
            if( i==1 ){
                localEnergy[b-1] += 0.25 * vertBonds[i-1];
            }
            else if( i==Lx-1 ){
                localEnergy[b-1] += 0.25 * vertBonds[i];
            }
        }
    } // for i

    return localEnergy;

}//localEnergy

std::vector<double> calculateLRenergy(int i, int Ly, MPS psi, std::vector<std::vector<ITensor>> LED_LR, SiteSet sites){

    std::vector<double> energy(Ly,0.);

    int index = (i-1)*Ly + 1;

    // long range interactions with smart ordering of gates 
    for(int m=0; m<=Ly-2; m++){

        psi.position(index+Ly+m);

        for(int n=Ly; n>m+1; n--){

            int b = index + n + m;
            auto g = BondGate(sites,b-1,b);
            auto AA = psi(b-1)*psi(b)*g.gate(); //contract over bond b
            AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
            psi.svdBond(g.i1(), AA, Fromright); //svd to restore MPS
            psi.position(g.i1()); //orthogonality center moves to the left

        } // for n

    } // for m
                
    for(int m = 0; m<Ly; m++){

        psi.position(index+2*m);
        auto ket = psi(index+2*m)*psi(index+2*m+1);
        energy[m] = eltC( dag(prime(ket,"Site")) * LED_LR[i-1][m] * ket).real();

    } // for m
    
    return energy;
}

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
double spinspin(int center, int b, MPS psi, SiteSet sites){
    
    double corrX;

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
        auto SxSx = 4.0*sites.op("Sx",center)*sites.op("Sx",center+1);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
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
        auto SxSx = 4.0*sites.op("Sx",center-1)*sites.op("Sx",center);
        corrX = eltC( dag(prime(ket,"Site")) * SxSx * ket).real();
    }
    else{
        corrX = 1.;
    }

    return corrX;

}//SxSx