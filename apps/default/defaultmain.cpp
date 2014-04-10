#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>

#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;

ElectronicSystem *setupSystem(string name);


int main(int argc, char **argv)
{

    int rank = 0;
#if USE_MPI
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    rank = world.rank();
#endif


    clock_t begin = clock();

    /********************************************************************************/

    //options:
    string method = "rhf";
    string chemicalSystem = "H2O";
    if(rank==0){

        cout << "---------------------------Hartree-Fock------------------------------"  << endl;
        cout << "system:    " << chemicalSystem << endl;
        cout << "method:    " << method << endl;
    }


    //Setup system:
    ElectronicSystem *system = setupSystem(chemicalSystem);


    //Choose method:
    HFsolver* solver;
    if(method == "rhf"){
        solver = new RHF(system);
    }else if(method == "uhf"){
        solver = new UHF(system);
    }else{
        cerr << "unknown method!" << endl;
        exit(0);
    }

    solver->runSolver();

    Analyser analyser(system,solver);
    analyser.calculateElectrostaticPotential();
//    analyser.atomicPartialCharge();
//    analyser.calculateChargeDensity();


    /********************************************************************************/
    clock_t end = clock();
    if(rank==0){
        cout << "Total elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;
    }

    return 0;

}

ElectronicSystem* setupSystem(string name)
{
    vector<Atom *> atoms;
    vector<rowvec3> atomsPos;

    if(name =="H2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -0.69, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 0.69, 0, 0 }));

    }else if(name =="Li2"){
        atomsPos.push_back({-2.5255, 0.0, 0.0});
        atomsPos.push_back({ 2.5255, 0.0, 0.0});
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-1.14, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 1.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        double D = 2.8;
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {D, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -D*cos((180-104.45) *M_PI/180.0),
                                                                              D*sin((180-104.45) *M_PI/180.0), 0.0}));
    }else if(name =="CO2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0, 0.0, 0.0}));

    }else if(name =="CH4"){        
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.043/sqrt(3), 2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-2.043/sqrt(3), -2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.043/sqrt(3), -2.043/sqrt(3), -2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-2.043/sqrt(3), 2.043/sqrt(3), -2.043/sqrt(3)}));

    }else if(name =="SiO4"){
        double D = 4.9;
        double T = 2.0*D /sqrt(3);
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {D/sqrt(3), D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), -D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{D/sqrt(3), -D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {T, -T, -T}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T - D/sqrt(3), -T -D/sqrt(3), -T - D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T + D/sqrt(3), -T +D/sqrt(3), -T -D/sqrt(3)}));

    }else if(name =="Fe2S2"){        
        double FeFe = 4.556129358;
        double SS   = 6.908838214;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0,  FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0, -FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {-SS*0.5, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {SS*0.5, 0.0, 0.0}));

    }else if(name =="FeS2"){
        double FeS = 3.802128689;
        double SFeS = 114.0;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31Gd.tm", {0.0, 0.0 , 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm", {FeS, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm",  {-FeS*cos((180.0-SFeS) *M_PI/180.0),
                                                                               FeS*sin((180.0-SFeS) *M_PI/180.0), 0.0}));

    }else if(name =="benzene"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        ,  4.68463073,  0.         }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,   0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023,  0.        }));


        double R = 6.666953288;
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.         ,  4.68463073,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,   R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023,  R        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023,  R        }));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    return system;

}

