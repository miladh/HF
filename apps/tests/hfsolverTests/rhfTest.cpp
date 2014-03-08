#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;
SUITE(DEVELOPMENT) {
TEST(H2_QZ)
{
    /*
     * test case:   H2
     * basis:       quadruple
     * bondlength:  1
     * energy:      -1.078547609
     *
     * source:
     *      Computational Physics
     *      Jos Thijssen
     * */

    int nElectrons;
    rowvec coreCharges,coreMass, A, B;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;

    nElectrons = 2;
    A = {-0.5, 0.0, 0.0};
    B = {0.5, 0.0, 0.0};
    coreCharges = {1 , 1};
    coreMass = {1 , 1};
    basisCoreA = new BasisSet("infiles/turbomole/H_Qzeta");
    basisCoreB = new BasisSet("infiles/turbomole/H_Qzeta");

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

    RHF *solver = new RHF(system,0,1);
    solver->runSolver();

    CHECK_CLOSE(-1.078547609, solver->getEnergy(), 1e-9);
}
TEST(H2_321G)
{
    int nElectrons;
    rowvec coreCharges,coreMass, A, B;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;

    nElectrons = 2;
    A = {-0.7, 0.0, 0.0};
    B = {0.7, 0.0, 0.0};
    coreCharges = {1 , 1};
    coreMass = {1 , 1};

    basisCoreA = new BasisSet("infiles/turbomole/H_3-21G");
    basisCoreB = new BasisSet("infiles/turbomole/H_3-21G");

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

    int nElectrons;
    rowvec coreCharges,coreMass, A, B;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;

    nElectrons = 2;
    A = {-0.69, 0.0, 0.0};
    B = {0.69, 0.0, 0.0};
    coreCharges = {1 , 1};
    coreMass = {1 , 1};

    basisCoreA = new BasisSet("infiles/turbomole/H_4-31G");
    basisCoreB = new BasisSet("infiles/turbomole/H_4-31G");

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

    int nElectrons;
    rowvec coreCharges,coreMass, A, B;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;

    nElectrons = 2;
    A = {-0.6925, 0.0, 0.0};
    B = {0.6925, 0.0, 0.0};
    coreCharges = {1 , 1};
    coreMass = {1 , 1};

    basisCoreA = new BasisSet("infiles/turbomole/H_6-31G_ds");
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

//    int nElectrons;
//    rowvec coreCharges,coreMass, A, B, C;

//    BasisSet *basisCoreA;
//    BasisSet *basisCoreB;
//    BasisSet *basisCoreC;


//    double R = 1.782;
//    double x = R * cos((180-104.45) *M_PI/180.0);
//    double y = R * sin((180-104.45) *M_PI/180.0);
//    //Water molecule
//    nElectrons = 10;
//    A = {R , 0.0, 0.0};
//    B = { -x, y, 0.0};
//    C = { 0.0, 0.0, 0.0};
//    coreCharges = {1 , 1, 8};
//    coreMass = {1 , 1, 16};

//    basisCoreA = new BasisSet("infiles/turbomole/H_6-31G_ds");
//    basisCoreB = new BasisSet("infiles/turbomole/H_6-31G_ds");
//    basisCoreC = new BasisSet("infiles/turbomole/O_6-31G_ds");

//    int maxAngularMomentum = basisCoreC->getAngularMomentum();
//    basisCoreA->setCoreCharge(coreCharges(0));
//    basisCoreA->setCoreMass(coreMass(0));
//    basisCoreA->setCorePosition(A);

//    basisCoreB->setCorePosition(B);
//    basisCoreB->setCoreCharge(coreCharges(1));
//    basisCoreB->setCoreMass(coreMass(1));

//    basisCoreC->setCorePosition(C);
//    basisCoreC->setCoreCharge(coreCharges(2));
//    basisCoreC->setCoreMass(coreMass(2));

//    System *system = new System(nElectrons, maxAngularMomentum);
//    system->addBasisSet(basisCoreA);
//    system->addBasisSet(basisCoreB);
//    system->addBasisSet(basisCoreC);

//    RHF *solver = new RHF(system,0,1);
//    solver->runSolver();

//    CHECK_CLOSE(-76.023551569545, solver->getEnergy(), 1e-9);

}
}
