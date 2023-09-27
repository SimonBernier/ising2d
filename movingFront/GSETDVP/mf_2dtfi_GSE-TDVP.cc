//
// Copyright [2023] [Simon Bernier]
//
#include "itensor/all.h"
#include "tdvp.h"
#include "basisextension.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, int, double, double, double, double, double);
//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<ITensor>,
                                        std::vector<std::vector<ITensor>>);
// calculate energy for periodic boundary at x = i
double calculatePBCenergy(int, int, MPS, ITensor, SiteSet);
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
    auto truncE = input.getReal("truncE", 1E-8);
    auto maxDim = input.getInt("maxDim", 128);
    auto tanhshift = input.getReal("tanhshift",2.0);
    auto dt = input.getReal("dt",0.1);

    // write results to file
    char schar[128];
    int n1 = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_v_%0.2f_tau_%0.1f_maxDim_%d_2dtfi_gsetdvp.dat",Ly,Lx,h,v,tau,maxDim);

    std::string s1(schar);
    std::ofstream datafile;
    datafile.open(s1); // opens the file
    if( !datafile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    datafile << "time" << " " << "en0(t)" << " " << "en-en0(t)" << " " << "svN" << " "
             << "localEn0" << " " << "localEn-localEn0" << " " << "sxsx0" << " " << "sxsx" << std::endl;

    auto N = Ly * Lx;
    auto sites = SpinHalf(N,{"ConserveQNs=",false,"ConserveParity",true});

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Lx, Ly, {"YPeriodic = ", true});

    // autompo hamiltonian
    for(auto j : lattice){
        ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    }
    // H at t=0 (gapped ground state)
    auto hvals = hvector(Lx, Ly, 0.0, h, v, tau, tanhshift);
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<=Ly; j++){
            int b = Ly*(i-1) + j;
            ampo += -2.0*hc*(1.0 + hvals[b-1]), "Sz", b;
        }
    }
    auto H = toMPO(ampo);

    //initial state
    auto state = InitState(sites); 
    for(auto j : range1(N)){
        state.set(j, "Up");
    }
    auto initState = MPS(state);

    // 2d ising model parameters
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 20, 50, 100, maxDim;
    sweeps.cutoff() = truncE;

    // calculate initial local energy density
    std::vector<double> localEn0(Ly*(Lx-1),0.0), localEn(Ly*(Lx-1),0.0); // local energy density vector
    std::vector<double> sxsx0(N,0.0), sxsx(N,0.0);

    //make 2D vector of ITensor for local energy operators
    //long-range interactions have the same structure as nearest-neighbour when we use swap gates
    std::vector<std::vector<ITensor>> LED(Lx, std::vector<ITensor>(Ly-1));
    std::vector<std::vector<ITensor>> LED_LR(Lx-1, std::vector<ITensor>(Ly));
    std::vector<std::vector<ITensor>> LED_OnSite(Lx, std::vector<ITensor>(Ly));
    std::vector<ITensor> LEDyPBC(Lx);
    // make local energy tensors
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<=Ly; j++){
            int index = (i-1)*Ly+j;
            //MPS long-range
            if(i<Lx && j==1){
                for(int m = 0; m<Ly; m++){
                    LED_LR[i-1][m] = -4.0*sites.op("Sx",index+2*m)*sites.op("Sx",index+2*m+1);                    
                }
            }
            // on-site operators
            LED_OnSite[i-1][j-1] = -2.0*hc*(1.0+hvals[index])*sites.op("Sz",index);
            //y-periodic boundary equations
            if(j==Ly){
                // site index-Ly+1 is moved to site index-1 with swap gates
                LEDyPBC[i-1] = -4.0*sites.op("Sx",index-1)*sites.op("Sx",index);
            }
            // MPS nearest-neighbour
            if(j<Ly){
                LED[i-1][j-1] = -4.0*sites.op("Sx",index)*sites.op("Sx",index+1);
            }
        }
    }

    //DMRG to find ground state at t=0
    auto [en0,psi0] = dmrg(H,initState,sweeps,{"Silent=",true});
    auto [en,psi] = dmrg(H,psi0,sweeps,{"Silent=",true});

    // calculate von Neumann S
    auto svN = vonNeumannS(psi, N/2);
    // calculate local energy density
    localEn = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LED_OnSite, LEDyPBC, LED_LR);
    for(int b = 1; b<=N; b++){
        sxsx[b-1] = spinspin(N/2+1, b, psi, sites);
    }

    // store data to file
    datafile << 0.0 << " " << en << " " << 0.0 << " " << svN << " ";
    for (int j=0; j < N-1; j++){
        datafile << localEn[j] << " ";
    }
    for (int j=0; j < N-1; j++){
        datafile << 0.0 << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << sxsx[j] << " ";
    }
    for (int j = 0; j < N; j++){
        datafile << sxsx[j] << " ";
    }
    datafile << std::endl;

    // time evolution parameters
    double tval = 0.0; //time
    double delta1 =  0.414490771794376*dt;
    double delta2 = -0.657963087177503*dt;
    double finalTime = 0.25*double(Lx) + 0.5*double(Lx)/v + 2.0*tau*tanhshift; // 0.5*N/2 + 0.5*N/v + 2*tau*shift
    int nt = int(finalTime/dt);

    // instantaneous GS search parameters
    sweeps = Sweeps(10); //number of sweeps is 5
    sweeps.maxdim() = maxDim; //gradually increase states kept
    sweeps.cutoff() = truncE; //desired truncation error
    // 4th order TDVP parameters
    auto sweeps1 = Sweeps(2); //two forward time steps of delta1
    sweeps1.maxdim() = maxDim;
    sweeps1.cutoff() = truncE;
    sweeps1.niter() = 20;
    auto sweeps2 = Sweeps(1); //one backward time step of delta2
    sweeps2.maxdim() = maxDim;
    sweeps2.cutoff() = truncE;
    sweeps2.niter() = 20;

    printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, en, svN, maxLinkDim(psi));

    ////////////////////////////////////////////////////////////////////////////
    ///////// time evolve //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for(int n=1; n<=nt; n++){

        tval += dt; //update time vector

	    // update magnetic field
        hvals = hvector(Lx, Ly, tval, h, v, tau, tanhshift);

        ampo = AutoMPO(sites);
        for(auto j : lattice){
            ampo += -4.0, "Sx", j.s1, "Sx", j.s2;
    	}
        for(int i=1; i<=Lx; i++){
            for(int j=1; j<=Ly; j++){
                int b = Ly*(i-1) + j;
                ampo += -2.0*hc*(1.0 + hvals[b-1]), "Sz", b;
            }
        }
        H = toMPO(ampo);

	    // calculate instantaneous ground state
        // use psi as initial condition
        en0 = dmrg(psi0,H,sweeps,{"Silent=",true});

        // time evolve with GSE-TDVP
        //std::vector<Real> epsilonK = {truncE, truncE};
        std::vector<int> dimK = {maxLinkDim(psi), maxLinkDim(psi)};
        //addBasis(psi, H, epsilonK, {"Cutoff",truncE,
        addBasis(psi, H, dimK, {"Cutoff",truncE,
                                        "Method", "DensityMatrix",
                                        "KrylovOrd",2,
                                        "Quiet",true});
        tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1});
        tdvp(psi, H, -Cplx_i*delta2, sweeps2, {"Silent",true,"Truncate",true,"NumCenter",1});
        auto en = tdvp(psi, H, -Cplx_i*delta1, sweeps1, {"Silent",true,"Truncate",true,"NumCenter",1});


        // change transverse fields
        for(int i=1; i<=Lx; i++){
            for(int j=1; j<=Ly; j++){
                int index = (i-1)*Ly+j;
                // on-site operators
                LED_OnSite[i-1][j-1] = -2.0*hc*(1.0+hvals[index])*sites.op("Sz",index);
            }
        }

        //calculate entanglement entropy
        svN = vonNeumannS(psi, N/2);
        // calculate local energy density <psi(t)|H(x,y)|psi(t)>
        localEn0 = calculateLocalEnergy(Lx, Ly, sites, psi0, LED, LED_OnSite, LEDyPBC, LED_LR);
        localEn = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LED_OnSite, LEDyPBC, LED_LR);
        for(int b = 1; b<=N; b++){
	        sxsx0[b-1] = spinspin(N/2+1, b, psi0, sites);
            sxsx[b-1] = spinspin(N/2+1, b, psi, sites);
        }

        datafile << tval << " " << en0 << " " << en-en0 << " " << svN << " ";
        for (int j = 0; j < N-1; j++){
            datafile << localEn0[j] << " ";
        }
        for (int j = 0; j < N-1; j++){
            datafile << localEn[j]-localEn0[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << sxsx0[j] << " ";
        }
        for (int j = 0; j < N; j++){
            datafile << sxsx[j] << " ";
        }
        datafile << std::endl;

        printfln("t = %0.2f, en-en0 = %0.3g, SvN = %0.3f, maxDim = %d", tval, en-en0, svN, maxLinkDim(psi));

    }

    datafile.close();

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
                                std::vector<std::vector<ITensor>> LED,
                                std::vector<std::vector<ITensor>> LED_OnSite,
                                std::vector<ITensor> LEDyPBC,
                                std::vector<std::vector<ITensor>> LED_LR){

    std::vector<std::vector<double>> tempEn(Lx, std::vector<double>(Ly, 0.0));

    std::vector<double> localEnergy(Ly*(Lx-1), 0.0); // interpolated energy density

    // MPS nearest-neighbour interaction
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<Ly; j++){ 

            int index = (i-1)*Ly + j;
            psi.position(index);
            auto ket = psi(index)*psi(index+1);
            tempEn[i-1][j-1] = eltC(dag(prime(ket,"Site")) * LED[i-1][j-1] * ket).real();

        }// for j
    }// for i

    // y-periodic interactions
    for(int i=1; i<=Lx; i++){

        tempEn[i-1][Ly-1] = calculatePBCenergy(i, Ly, psi, LEDyPBC[i-1], sites);

    }// for i

    // MPS long-range interactions
    for(int i=1; i<Lx; i++){

        int index = (i-1)*Ly + 1;
        auto lrEnergy = calculateLRenergy(i, Ly, psi, LED_LR, sites);

        for(int m = 0; m<Ly; m++){ // long-range

            localEnergy[index-1+m] = lrEnergy[m];

        }// for m
    }// for i

    // MPS on-site interactions
    for(int i=1; i<Lx; i++){
        for(int j = 1; j <= Ly; j++){

            int index = (i-1)*Ly + j;

            psi.position(index);
            if(i == 1){
                localEnergy[index-1] += eltC(dag(prime(psi(index),"Site")) * LED_OnSite[i-1][j-1] * psi(index)).real();
            }
            else{
                localEnergy[index-1] += 0.5*eltC(dag(prime(psi(index),"Site")) * LED_OnSite[i-1][j-1] * psi(index)).real();
            }
            
            psi.position(index+Ly);
            if(i==Lx-1){
                localEnergy[index-1] += eltC(dag(prime(psi(index+Ly),"Site")) * LED_OnSite[i][j-1] * psi(index+Ly)).real();
            }
            else{
                localEnergy[index-1] += 0.5*eltC(dag(prime(psi(index+Ly),"Site")) * LED_OnSite[i][j-1] * psi(index+Ly)).real();
            }
        }// for j
    }// for i

    // interpolate tempEn into localEnergy
    for(int i=1; i<Lx; i++){
        for(int j=1; j<=Ly; j++){ 

            int index = (i-1)*Ly + j;

            if(j == 1){ // interpolate with periodic boundary conditions
                localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-1] + tempEn[i][j-1] + tempEn[i-1][j+Ly-2] + tempEn[i][j+Ly-2] );
            }
            else{ // interpolate normally
                localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-2] + tempEn[i][j-2] + tempEn[i-1][j-1] + tempEn[i][j-1] );
            }
            // add left/right boundary terms
            if( i==1 ){
                if( j== 1){
                    localEnergy[index-1] += 0.25 * ( tempEn[i-1][j+Ly-2] + tempEn[i-1][j-1]);
                }
                else{
                    localEnergy[index-1] += 0.25 * ( tempEn[i-1][j-2] + tempEn[i-1][j-1]);
                }
            }// if i=1
            else if( i==Lx-1 ){
                if( j== 1){
                    localEnergy[index-1] += 0.25 * ( tempEn[i][j+Ly-2] + tempEn[i][j-1]);
                }
                else{
                    localEnergy[index-1] += 0.25 * ( tempEn[i][j-2] + tempEn[i][j-1]);
                }
            } // if i=Lx-1
        } // for j
    } // for i

    return localEnergy;

}//localEnergy

double calculatePBCenergy(int i, int Ly, MPS psi, ITensor hterm, SiteSet sites){

    double energy;
    int index = (i-1)*Ly + 1;

    psi.position(index);

    for(int n=0; n<Ly-2; n++){

        int b = index+n;//define gate bond
        auto g = BondGate(sites,b,b+1);
        auto AA = psi(b)*psi(b+1)*g.gate(); //contract over bond b
        AA.replaceTags("Site,1","Site,0"); //replace site tags for correct svd
        psi.svdBond(g.i1(), AA, Fromleft); //svd to restore MPS
        psi.position(g.i2()); //orthogonality center moves to the right

    } // for n

    auto ket = psi(index+Ly-2)*psi(index+Ly-1);
    energy = eltC( dag(prime(ket,"Site")) * hterm * ket).real();

    return energy;
} // calculatePBCenergy

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
