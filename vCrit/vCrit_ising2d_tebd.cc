#include "itensor/all.h"

using namespace itensor;

//makes gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int, int, double, SiteSet, std::vector<std::vector<ITensor>>,
                                std::vector<ITensor>, std::vector<std::vector<ITensor>>);
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
    auto h = input.getReal("h", 3.0);
    auto truncE = input.getReal("truncE", 1E-10);
    auto maxB = input.getInt("maxB", 128);
    
    printfln("Ly = %d, Lx = %d, h = %0.2f", Ly, Lx, h);

    // write results to file
    char schar[64];
    int n = std::sprintf(schar,"Ly_%d_Lx_%d_h_%0.2f_vCrit.dat",Ly,Lx,h);
    std::string s(schar);
    std::ofstream dataFile;
    dataFile.open(s); // opens the file
    if( !dataFile ) { // file couldn't be opened
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }
    //make header for t=0 calculations
    dataFile << "t=0" << " " << "enPsi" << " " << "MaxDimPsi" << " " << "enPhi" << " " << "MaxDimPhi" << " " 
             << "SvN(x,t=0)" << " " << std::endl;

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

    //initial state
    auto initState = InitState(sites); 
    for(auto j : range1(N)){
        initState.set(j, "Up");
    }

    // 2d ising model parameters
    auto sweeps = Sweeps(15);
    sweeps.maxdim() = 10, 20, 100, 100, 200, 200, 400, 400, 512;
    sweeps.cutoff() = 1E-10;
    sweeps.noise() = 1E-7,1E-8,0.0;
    
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
                LEDyPBC[i-1] =  -4.0  *sites.op("Sx",index-1)*sites.op("Sx",index);
                LEDyPBC[i-1] += -2.0*h*sites.op("Id",index-1)*sites.op("Sz",index);
            }
            // MPS nearest-neighbour
            if(j<Ly){
                LED[i-1][j-1] =  -4.0  *sites.op("Sx",index)*sites.op("Sx",index+1);
                LED[i-1][j-1] += -2.0*h*sites.op("Sz",index)*sites.op("Id",index+1);
            }
        }
    }
        
    // calculate ground state
    auto H = toMPO(ampo);
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

    printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",0,0., en_phi,maxLinkDim(phi));
    // store to file
    dataFile << 0.0 << " " << en_psi << " " << maxLinkDim(psi) << " " << en_phi << " " << maxLinkDim(phi) << " ";
    for(int j = 0; j<Lx-1; j++){ //save svn
        dataFile << svn[j] << " ";
    }
    dataFile << std::endl;

    dataFile << "tval" << " " << "enPhi" << " " << "MaxDimPhi" << " " << "SvN(x,t)" << " " << std::endl;

    // time evolution parameters
    double tval = 0.0; //time
    double ttotal = 2.0;
    double tstep = 0.1;
    Real dt = 0.01; //time step
    int Nt = int(ttotal/tstep); //number of time steps

    printfln("\nStarting second order TEBD, dt = %0.2f\n", dt);
    Args args = Args("Cutoff=",1E-10,"MaxDim=",512);

    std::vector<BondGate> gates = makeGates(Lx, Ly, dt, sites, LED, LEDyPBC, LED_LR);

    ////////////////////////////////////////////////////////////////////////////
    ///////// time evolve //////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    for(int n=1; n<=Nt; n++){
        tval += tstep; //update time vector
        
        //4th order TEBD
        gateTEvol(gates,tstep,dt,phi,{args,"Verbose=",false});
        phi.orthogonalize(args);

        en_phi = innerC(phi, H, phi).real();

        for(auto i : range1(Lx-1)){
            auto b = i*Ly;
            svn[i-1] = vonNeumannS(phi, b);
        }

        //write to file
        dataFile << tval << " " << en_phi << " " << maxLinkDim(phi) << " ";
        for(int j = 0; j<Lx-1; j++){ //save svn
            dataFile << svn[j] << " ";
        }
        dataFile << std::endl;

        printfln("\nIteration %d, time = %0.2f; phi energy = %0.3f, max link dim is %d",n,tval,en_phi,maxLinkDim(phi));

    }// for n

    dataFile.close();

    print(" END PROGRAM. TIME TAKEN :");
    printfln("Time taken: %.3fs\n", (double)(std::clock() - tStart)/CLOCKS_PER_SEC);

    return 0;

}//main

// second order Trotter breakup of time step dt
// returns a vector of gates to pass to function gateTEvol
std::vector<BondGate> makeGates(int Lx, int Ly, double dt, SiteSet sites, 
                                std::vector<std::vector<ITensor>> LED,
                                std::vector<ITensor> LEDyPBC,
                                std::vector<std::vector<ITensor>> LED_LR){

    std::vector<BondGate> gates; 
    
    //Create the gates

    // vertical bonds
    for(int i=1; i<=Lx; i++){
        for(int j=1; j<=Ly; j++){ 
            int index = (i-1)*Ly + j; //MPS site index

            // nearest neighbour
            if(j<Ly){
                auto hterm = LED[i-1][j-1];
                auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
                gates.push_back(g);
            }// if j<Ly
            
            // y-PBC
            if(j==Ly){
                int index = (i-1)*Ly + 1; //MPS site index
                for(int n=0; n<Ly-2; n++){ //swap from index-Ly+1 to index-1
                    int b = index+n;
                    auto swapGate = BondGate(sites,b,b+1);
                    gates.push_back(swapGate);
                }// for n

                auto hterm = LEDyPBC[i-1];
                auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
                gates.push_back(g);

                //restore the state to the original MPS
                for(int n=Ly-2; n>0; n--){
                    int b = index+n;
                    auto swapGate = BondGate(sites,b-1,b);
                    gates.push_back(swapGate);
                }// for n
            } // if j==1
        } // for j
    } // for i

    // long-range interaction
    for(int i=1; i<Lx; i++){

        int index = (i-1)*Ly + 1; //MPS site index

        for(int m=0; m<=Ly-2; m++){
            for(int n=Ly; n>m+1; n--){
                int b = index + n + m;
                auto swapGate = BondGate(sites,b-1,b);
                gates.push_back(swapGate);
            }
        }

        for(int m = 0; m<Ly; m++){
            auto hterm = LED_LR[i-1][m];
            auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }

        // bring index+1 back to position index+Ly
        for(int m=Ly-2; m>=0; m--){
            for(int n=m+1; n<Ly; n++){
                int b = index + n + m;
                auto swapGate = BondGate(sites,b,b+1);
                gates.push_back(swapGate);
            }
        }
    }// for i

    // make reverse gates

    // reverse long-range
    for(int i=Lx-1; i>=1; i--){

        int index = (i-1)*Ly + 1; //MPS site index

        for(int m=0; m<=Ly-2; m++){
            for(int n=Ly; n>1+m; n--){
                int b = index + n + m;
                auto swapGate = BondGate(sites,b-1,b);
                gates.push_back(swapGate);
            }
        }
              
        for(int m = Ly-1; m>=0; m--){
            auto hterm = LED_LR[i-1][m];
            auto g = BondGate(sites,index+2*m,index+2*m+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
        }

        // bring index+1 back to position index+Ly
        for(int m=Ly-2; m>=0; m--){
            for(int n=m+1; n<Ly; n++){
                int b = index + n + m;
                auto swapGate = BondGate(sites,b,b+1);
                gates.push_back(swapGate);
            }
        }
    }// for i

    // reverse vertical bonds
    for(int i=Lx; i>=1; i--){
        for(int j=Ly; j>=1; j--){ 

            int index = (i-1)*Ly + j; //MPS site index

            // nearest-neighbour
            if(j<Ly){
            auto hterm = LED[i-1][j-1];
            auto g = BondGate(sites,index,index+1,BondGate::tReal,dt/2.,hterm);
            gates.push_back(g);
            } // if j<Ly

            // y-PBC
            if(j==Ly){
                int index = (i-1)*Ly + 1; //MPS site index
                for(int n=0; n<Ly-2; n++){ //swap from index-Ly+1 to index-1
                    int b = index+n;
                    auto swapGate = BondGate(sites,b,b+1);
                    gates.push_back(swapGate);
                }// for n
                    
                auto hterm = LEDyPBC[i-1];
                auto g = BondGate(sites,index+Ly-2,index+Ly-1,BondGate::tReal,dt/2.,hterm);
                gates.push_back(g);

                //restore the state to the original MPS
                for(int n=Ly-2; n>0; n--){
                    int b = index+n;
                    auto swapGate = BondGate(sites,b-1,b);
                    gates.push_back(swapGate);
                }// for n
            } // if j==1
        } // for j
    } // for i

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