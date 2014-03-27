#include <unittest++/UnitTest++.h>
#include <hf.h>
#include <mpi.h>
#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;
SUITE(DEVELOPMENT) {
    TEST(H2_321G)
    {
        /*
     * test case:   H2
     * basis:       3-21G
     * bondlength:  1
     * energy:      -1.122933364
     *
     * source:
     *      unknown
     *
     * */


        if(MPI::COMM_WORLD.Get_rank()==0){
            cout << "system:    " << "H2" << endl;
            cout << "method:    " << "RHF" << endl;
            cout << "basis:     " << "3-21G" << endl;
        }

        Atom *atomA;
        Atom *atomB;

        atomA = new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-0.7, 0.0, 0.0});
        atomB = new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 0.7, 0.0, 0.0});
        int maxAngularMomentum = atomA->angularMomentum();

        ElectronicSystem *system = new ElectronicSystem (maxAngularMomentum);
        system->addAtom(atomA);
        system->addAtom(atomB);

        RHF *solver = new RHF(system,0,1);
        solver->runSolver();

        CHECK_CLOSE(-1.122933364, solver->getEnergy(), 1e-9);
    }

    TEST(H2_431G)
    {
        /*
     * test case:   H2
     * basis:       4-31G
     * bondlength:  1.380
     * energy:      -1.127
     *
     * source:
     *      Molecular Quantum Mechanics
     *      Peter Atkins
     * */

        if(MPI::COMM_WORLD.Get_rank()==0){
            cout << "system:    " << "H2" << endl;
            cout << "method:    " << "RHF" << endl;
            cout << "basis:     " << "4-31G" << endl;
        }


        Atom *atomA;
        Atom *atomB;

        atomA = new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", {-0.69, 0.0, 0.0});
        atomB = new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { 0.69, 0.0, 0.0});
        int maxAngularMomentum = atomA->angularMomentum();

        ElectronicSystem *system = new ElectronicSystem (maxAngularMomentum);
        system->addAtom(atomA);
        system->addAtom(atomB);

        RHF *solver = new RHF(system,0,1);
        solver->runSolver();

        CHECK_CLOSE(-1.12682776, solver->getEnergy(), 1e-9);
    }
    TEST(H2_631G_ds)
    {
        /*
     * test case:   H2
     * basis:       6-31G**
     * bondlength:  1.385
     * energy:      -1.131
     *
     * source:
     *      Molecular Quantum Mechanics
     *      Peter Atkins
     * */

        if(MPI::COMM_WORLD.Get_rank()==0){
            cout << "system:    " << "H2" << endl;
            cout << "method:    " << "RHF" << endl;
            cout << "basis:     " << "6-31G**" << endl;
        }

        Atom *atomA;
        Atom *atomB;

        atomA = new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {-0.6925, 0.0, 0.0});
        atomB = new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { 0.6925, 0.0, 0.0});
        int maxAngularMomentum = atomA->angularMomentum() + 1;

        ElectronicSystem *system = new ElectronicSystem (maxAngularMomentum);
        system->addAtom(atomA);
        system->addAtom(atomB);

        RHF *solver = new RHF(system,0,1);
        solver->runSolver();

        CHECK_CLOSE(-1.1313335068087, solver->getEnergy(), 1e-9);

    }

    TEST(H2O_431G)
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

        if(MPI::COMM_WORLD.Get_rank()==0){
            cout << "system:    " << "H2O" << endl;
            cout << "method:    " << "RHF" << endl;
            cout << "basis:     " << "4-31G" << endl;
        }


        Atom *atomA;
        Atom *atomB;
        Atom *atomC;

        double x = 1.797*cos((180-104.45) *M_PI/180.0);
        double y = 1.797*sin((180-104.45) *M_PI/180.0);

        atomA = new Atom("infiles/turbomole/atom_8_basis_4-31G.tm", { 0.0, 0.0, 0.0});
        atomB = new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", {1.797, 0.0, 0.0});
        atomC = new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -x, y, 0.0});
        int maxAngularMomentum = atomA->angularMomentum();

        ElectronicSystem *system = new ElectronicSystem (maxAngularMomentum);
        system->addAtom(atomA);
        system->addAtom(atomB);
        system->addAtom(atomC);

        RHF *solver = new RHF(system,0,1);
        solver->runSolver();

        CHECK_CLOSE(-75.907340813845, solver->getEnergy(), 1e-9);

    }
}

SUITE(SLOWTESTS) {
    TEST(H2O_631Gds)
    {
        /*
     * test case:   H2O
     * basis:       6-31G**
     * bondlength:  1.782
     * bond angle:  104.45
     * energy:      -76.023
     *
     * source:
     *      Molecular Quantum Mechanics
     *      Peter Atkins
     * */

        if(MPI::COMM_WORLD.Get_rank()==0){
            cout << "system:    " << "H2O" << endl;
            cout << "method:    " << "RHF" << endl;
            cout << "basis:     " << "6-31G**" << endl;
        }


        Atom *atomA;
        Atom *atomB;
        Atom *atomC;

        double R = 1.782;
        double x = R * cos((180-104.45) *M_PI/180.0);
        double y = R * sin((180-104.45) *M_PI/180.0);

        atomA = new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", { 0.0, 0.0, 0.0});
        atomB = new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", { -x, y, 0.0});
        atomC = new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {R , 0.0, 0.0});
        int maxAngularMomentum = atomA->angularMomentum();

        ElectronicSystem *system = new ElectronicSystem (maxAngularMomentum);
        system->addAtom(atomA);
        system->addAtom(atomB);
        system->addAtom(atomC);

        RHF *solver = new RHF(system,0,1);
        solver->runSolver();

        CHECK_CLOSE(-76.023551569545, solver->getEnergy(), 1e-9);

    }
}
