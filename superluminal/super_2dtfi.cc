#include "itensor/all.h"
#include "tdvp.h"

using namespace itensor;

//magnetic field vector
std::vector<double> hvector(int, int, double, double, double, double, double);
//local energy calculation using swap gates
std::vector<double> calculateLocalEnergy(int, int, SiteSet, MPS, 
                                        std::vector<std::vector<ITensor>>,
                                        std::vector<ITensor>,  
                                        std::vector<std::vector<ITensor>>);
// calculate energy for periodic boundary at x = i
double calculatePBCenergy(int, int, MPS, ITensor, SiteSet);
// calculate energy of horizontal bonds for one column
std::vector<double> calculateLRenergy(int, int, MPS, std::vector<std::vector<ITensor>>, SiteSet);
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
    std::ofstream enerfile;
    enerfile.open(s1); // opens the file
    if( !enerfile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header
    enerfile << "tval" << " " << "energy" << " " << "svn" << " " << "MaxDim" << " " << "localEnergy" << " " << std::endl;

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
  
    //make 2D vector of ITensor for local energy operators
    //long-range interactions have the same structure as nearest-neighbour when we use swap gates
    std::vector<std::vector<ITensor>> LED(Lx, std::vector<ITensor>(Ly-1));
    std::vector<std::vector<ITensor>> LED_LR(Lx-1, std::vector<ITensor>(Ly));
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
            //y-periodic boundary equations
            if(j==Ly){
                // site index-Ly+1 is moved to site index-1 with swap gates
                LEDyPBC[i-1] = -4.0*sites.op("Sx",index-1)*sites.op("Sx",index);
                LEDyPBC[i-1] += -hc*sites.op("Id",index-1)*sites.op("Sz",index);
                LEDyPBC[i-1] += -hc*sites.op("Sz",index-1)*sites.op("Id",index);
            }
            // MPS nearest-neighbour
            if(j<Ly){
                LED[i-1][j-1] = -4.0*sites.op("Sx",index)*sites.op("Sx",index+1);
                LED[i-1][j-1] += -hc*sites.op("Sz",index)*sites.op("Id",index+1);
                LED[i-1][j-1] += -hc*sites.op("Id",index)*sites.op("Sz",index+1);
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
    localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);

    //store to file
    enerfile << 0.0 << " " << energy << " " << svN << " " << maxLinkDim(psi) << " "; //print to file
    for(int j = 0; j<Ly*(Lx-1); j++){ //save local energy values
        enerfile << localEnergy[j] << " ";
    }
    enerfile << std::endl;

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
        localEnergy = calculateLocalEnergy(Lx, Ly, sites, psi, LED, LEDyPBC, LED_LR);

        //write to file
        enerfile << tval << " " << energy << " " << svN << " " << maxLinkDim(psi) << " ";
        for(int j = 0; j<Ly*(Lx-1); j++){ //save local energy values
            enerfile << localEnergy[j] << " ";
        }
        enerfile << std::endl;

        printfln("t = %0.2f, energy = %0.3f, SvN = %0.3f, maxDim = %d", tval, energy, svN, maxLinkDim(psi));

    }

    enerfile.close();

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

        for(int m = 0; m<Ly; m++){

            localEnergy[index-1+m] = lrEnergy[m];

        } // for m

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

    //smart ordering of gates 
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
