#include <unittest++/UnitTest++.h>
#include <hf.h>
#include<mpi.h>
#include <armadillo>
#include <iostream>
#include <fstream>


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

    if(MPI::COMM_WORLD.Get_rank()==0){
        cout << "system:    " << "H2O" << endl;
        cout << "method:    " << "UHF" << endl;
        cout << "basis:     " << "4-31G" << endl;
    }

    int nElectrons;
    rowvec coreCharges,coreMass, A, B, C;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;
    BasisSet *basisCoreC;

    double x = 1.797*cos((180-104.45) *M_PI/180.0);
    double y = 1.797*sin((180-104.45) *M_PI/180.0);
    //Water molecule
    nElectrons = 10;
    A = {1.797, 0.0, 0.0};
    B = { -x, y, 0.0};
    C = { 0.0, 0.0, 0.0};
    coreCharges = {1 , 1, 8};
    coreMass = {1 , 1, 16};

    basisCoreA = new BasisSet("infiles/turbomole/H_4-31G");
    basisCoreB = new BasisSet("infiles/turbomole/H_4-31G");
    basisCoreC = new BasisSet("infiles/turbomole/O_4-31G");

    int maxAngularMomentum = basisCoreC->getAngularMomentum();
    basisCoreA->setCoreCharge(coreCharges(0));
    basisCoreA->setCoreMass(coreMass(0));
    basisCoreA->setCorePosition(A);

    basisCoreB->setCorePosition(B);
    basisCoreB->setCoreCharge(coreCharges(1));
    basisCoreB->setCoreMass(coreMass(1));

    basisCoreC->setCorePosition(C);
    basisCoreC->setCoreCharge(coreCharges(2));
    basisCoreC->setCoreMass(coreMass(2));

    System *system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);
    system->addBasisSet(basisCoreC);

    UHF *solver = new UHF(system,0,1);
    solver->runSolver();

    CHECK_CLOSE(-75.907340813846, solver->getEnergy(), 1e-9);
}
}


SUITE(SLOWTESTS) {
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

    if(MPI::COMM_WORLD.Get_rank()==0){
        cout << "system:    " << "HF" << endl;
        cout << "method:    " << "UHF" << endl;
        cout << "basis:     " << "6-31G**" << endl;
    }

    int nElectrons;
    rowvec coreCharges,coreMass, A, B;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;

    nElectrons = 10;
    A = { -1.889725989, 0.0, 0.0};
    B = { 1.889725989, 0.0, 0.0};
    coreCharges = {9 , 1};
    coreMass = {19 , 1};

    basisCoreA = new BasisSet("infiles/turbomole/F_6-31G_ds");
    basisCoreB = new BasisSet("infiles/turbomole/H_6-31G_ds");

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCoreCharge(coreCharges(0));
    basisCoreA->setCoreMass(coreMass(0));
    basisCoreA->setCorePosition(A);

    basisCoreB->setCorePosition(B);
    basisCoreB->setCoreCharge(coreCharges(1));
    basisCoreB->setCoreMass(coreMass(1));

    System *system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);


    UHF *solver = new UHF(system,0,1);
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

    if(MPI::COMM_WORLD.Get_rank()==0){
        cout << "system:    " << "CH4" << endl;
        cout << "method:    " << "UHF" << endl;
        cout << "basis:     " << "6-31G*" << endl;
    }

    int nElectrons;
    rowvec coreCharges,coreMass, A, B,C,D,E;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;
    BasisSet *basisCoreC;
    BasisSet *basisCoreD;
    BasisSet *basisCoreE;

    nElectrons = 10;
    A = {0.0, 0.0, 0.0};
    B = {1.0, 1.0, 1.0};
    C = {-1.0, -1.0, 1.0};
    D = {1.0, -1.0, -1.0};
    E = {-1.0, 1.0, -1.0};

    double d =1.889725989;
    B *=2.052242424/sqrt(3);C *=2.052242424/sqrt(3);D *=2.052242424/sqrt(3);E *=d/sqrt(3);

    coreCharges = {6, 1 , 1, 1, 1};
    coreMass = {6 , 1, 1, 1, 1};


    basisCoreA = new BasisSet("infiles/turbomole/C_6-31G_d");
    basisCoreB = new BasisSet("infiles/turbomole/H_6-31G");
    basisCoreC = new BasisSet("infiles/turbomole/H_6-31G");
    basisCoreD = new BasisSet("infiles/turbomole/H_6-31G");
    basisCoreE = new BasisSet("infiles/turbomole/H_6-31G");

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCoreCharge(coreCharges(0));
    basisCoreA->setCoreMass(coreMass(0));
    basisCoreA->setCorePosition(A);

    basisCoreB->setCorePosition(B);
    basisCoreB->setCoreCharge(coreCharges(1));
    basisCoreB->setCoreMass(coreMass(1));

    basisCoreC->setCorePosition(C);
    basisCoreC->setCoreCharge(coreCharges(2));
    basisCoreC->setCoreMass(coreMass(2));

    basisCoreD->setCorePosition(D);
    basisCoreD->setCoreCharge(coreCharges(3));
    basisCoreD->setCoreMass(coreMass(3));

    basisCoreE->setCorePosition(E);
    basisCoreE->setCoreCharge(coreCharges(4));
    basisCoreE->setCoreMass(coreMass(4));

    System *system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);
    system->addBasisSet(basisCoreC);
    system->addBasisSet(basisCoreD);
    system->addBasisSet(basisCoreE);

    UHF *solver = new UHF(system,0,1);
    solver->runSolver();

    double energy = -40.349369;
    double dE = abs(solver->getEnergy() - energy);
    CHECK_CLOSE(0.159802, dE, 1e-6);
}
}
