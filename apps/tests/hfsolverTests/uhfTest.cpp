#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <boost/mpi.hpp>
#include <hf.h>


using namespace std;
using namespace arma;
using namespace hf;

SUITE(DEVELOPMENT) {
    TEST(H20_431G)
    {
        /*
     * test case:   H2O
     * basis:       4-31G
     * bondlength:  1.797
     * bond angle:  104.45
     * energy:      -75.907
     *
     * source:
     *      Molecular Quantum Mechanics
     *      Peter Atkins
     * */

        int myRank = 0;
#if USE_MPI
        boost::mpi::communicator world;
        myRank = world.rank();
#endif
        if(myRank==0){
            cout << "system:    " << "H2O" << endl;
            cout << "method:    " << "UHF" << endl;
            cout << "basis:     " << "4-31G" << endl;
        }

        vector<Atom *> atoms;
        double x = 1.797*cos((180-104.45) *M_PI/180.0);
        double y = 1.797*sin((180-104.45) *M_PI/180.0);

        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_4-31G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", {1.797, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -x, y, 0.0}));

        ElectronicSystem *system = new ElectronicSystem ();
        system->addAtoms(atoms);

        UHF *solver = new UHF(system);
        solver->runSolver();

        CHECK_CLOSE(-75.907340813846, solver->getEnergy(), 1e-9);
    }
}


SUITE(SLOWTESTS_UHF) {
    TEST(HF_6_31G_ds)
    {
        /*
     * test case:    HF
     * basis:        6-31G**
     * bondlength:   3.778
     * energy(FCI): -100.020326
     * Error (UHF):  0.154528
     *
     * source:
     *      J. Chem. Phys. 118, 1610 (2003);
     *      http://dx.doi.org/10.1063/1.1531658
     * */
        int myRank = 0;
#if USE_MPI
        boost::mpi::communicator world;
        myRank = world.rank();
#endif
        if(myRank==0){
            cout << "system:    " << "HF" << endl;
            cout << "method:    " << "UHF" << endl;
            cout << "basis:     " << "6-31G**" << endl;
        }

        vector<Atom *> atoms;
        atoms.push_back( new Atom("infiles/turbomole/atom_9_basis_6-31Gds.tm", { -1.889725989, 0.0, 0.0}));
        atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { 1.889725989, 0.0, 0.0}));

        ElectronicSystem *system = new ElectronicSystem ();
        system->addAtoms(atoms);

        UHF *solver = new UHF(system);
        solver->runSolver();

        double E = -100.020326;
        double dE = abs(solver->getEnergy() - E);
        CHECK_CLOSE(0.154528, dE, 1e-6);

    }


    TEST(CH4_6_31G_d)
    {
        /*
     * test case:    CH4
     * basis:        6-31G*
     * bondlength:   1.889725989 (2.052242424)
     * energy(FCI): -40.349369
     * Error (UHF):  0.159802
     *
     * source:
     *      J. Chem. Phys. 118, 1610 (2003);
     *      http://dx.doi.org/10.1063/1.1531658
     * */

        int myRank = 0;
#if USE_MPI
        boost::mpi::communicator world;
        myRank = world.rank();
#endif
        if(myRank==0){
            cout << "system:    " << "CH4" << endl;
            cout << "method:    " << "UHF" << endl;
            cout << "basis:     " << "6-31G*" << endl;
        }

        rowvec A, B,C,D,E;
        A = {0.0, 0.0, 0.0};
        B = {1.0, 1.0, 1.0};
        C = {-1.0, -1.0, 1.0};
        D = {1.0, -1.0, -1.0};
        E = {-1.0, 1.0, -1.0};

        double d =1.889725989;
        B *=2.052242424/sqrt(3);C *=2.052242424/sqrt(3);D *=2.052242424/sqrt(3);E *=d/sqrt(3);

        vector<Atom *> atoms;
       atoms.push_back( new Atom("infiles/turbomole/atom_6_basis_6-31Gd.tm", A));
       atoms.push_back(  new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", B));
       atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", C));
       atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", D));
       atoms.push_back( new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", E));

        ElectronicSystem *system = new ElectronicSystem ();
        system->addAtoms(atoms);
        UHF *solver = new UHF(system);

        solver->runSolver();

        double energy = -40.349369;
        double dE = abs(solver->getEnergy() - energy);
        CHECK_CLOSE(0.159802, dE, 1e-6);
    }
}
