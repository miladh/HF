#include <unittest++/UnitTest++.h>
#include<system/system.h>
#include<hfSolver/hfsolver.h>

#include <armadillo>
#include <iostream>
#include <fstream>


using namespace std;
using namespace arma;

TEST(H2_QZ)
{
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

    HFsolver solver(system);
    solver.runSolver();

    CHECK_CLOSE(-1.078547609, solver.getEnergy(), 1e-9);
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

    HFsolver solver(system);
    solver.runSolver();

    CHECK_CLOSE(-1.122933364, solver.getEnergy(), 1e-9);

}
TEST(H2_431G)
{
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

    HFsolver solver(system);
    solver.runSolver();

    CHECK_CLOSE(-1.12682776, solver.getEnergy(), 1e-9);

}
TEST(H2O_431G)
{
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

    HFsolver solver(system);
    solver.runSolver();

    CHECK_CLOSE(-75.907340813845, solver.getEnergy(), 1e-9);

}
