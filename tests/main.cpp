#include <unittest++/UnitTest++.h>
#include <integrator/integrator.h>
#include<bomd/bomd.h>
#include<basisSet/splitValence/h_321g.h>

TEST(Hmolecules)
{
    BasisSet *basisCoreA, *basisCoreB, *basisCoreC, *basisCoreD;
    rowvec coreCharges,coreMass;
    rowvec A, B, C, D;

    int nElectrons = 4;
    A = {-0.7, 0.0, 0.0};
    B = {0.7, 0.0, 0.0};
    C = {0.0, 0.7, 2.0};
    D = {0.0, -0.7, -2.0};
    coreCharges = {1 , 1, 1, 1};
    coreMass = {1 , 1, 1, 1};

    basisCoreA  = new H_321G;
    basisCoreB  = new H_321G;
    basisCoreC  = new H_321G;
    basisCoreD  = new H_321G;


    int maxAngularMomentum = basisCoreA->getAngularMomentum()+1;
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

    System *system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);
    system->addBasisSet(basisCoreC);
    system->addBasisSet(basisCoreD);

    BOMD boSolver(system);
    boSolver.runDynamics();

}

int main()
{
    return UnitTest::RunAllTests();
}
