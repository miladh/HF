#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>

#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;

ElectronicSystem* setupSystem(string name);
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
        cout << "---------------------------BOMD------------------------------"  << endl;
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


    BOMD boSolver(system, solver);
    boSolver.runDynamics();

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
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 1, 0, 0 }));

    }else if(name =="Li2"){
        atomsPos.push_back({-2.5255, 0.0, 0.0});
        atomsPos.push_back({ 2.5255, 0.0, 0.0});
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-1.14, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 1.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {1.797, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1.797*cos((180-104.45) *M_PI/180.0),
                                                                              1.797*sin((180-104.45) *M_PI/180.0), 0.0}));
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
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_3-21G.tm", {0.0, 1.0 , 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_3-21G.tm", {1.0 , 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21G.tm", {1.0, 1.0, 0.0}));

    }else if(name =="benzene"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.99261000, 0.99261000, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-1.35593048, 0.36332048, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {0.36332048, -1.35593048, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-0.99261000, -0.99261000, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {1.35593048, -0.36332048, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-0.36332048, 1.35593048, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {1.75792000, 1.75792000, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-2.40136338, 0.64344338, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.64344338, -2.40136338, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-1.75792000, -1.75792000, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {2.40136338, -0.64344338, 0.00000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-0.64344338, 2.40136338, 0.00000000}));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

    return system;


}

