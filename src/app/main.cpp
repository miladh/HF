#include <iostream>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<contractedGTO/contractedGTO.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/h_321g.h>
#include<basisSet/o_321g.h>


using namespace arma;
using namespace std;


int main()
{

//    H_QuadZeta *basisCoreA  = new H_QuadZeta;
//    H_QuadZeta *basisCoreB  = new H_QuadZeta;

//    H_321G *basisCoreA = new H_321G;
//    H_321G *basisCoreB = new H_321G;

    O_321G *basisCoreA = new O_321G;
    O_321G *basisCoreB = new O_321G;

    int nBasisFunc = 9;
    int nNuclei    = 2;
    rowvec coreCharges = {8.0 , 8.0};

    rowvec A = {-1.41 , 0, 0};
    rowvec B = {1.41 , 0, 0};

    /********************************************************/



    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    int nOrbitals = nBasisFunc * nNuclei;

    basisCoreA->setCorePosition(A);
    basisCoreB->setCorePosition(B);


    System system(nOrbitals,maxAngularMomentum, coreCharges);
    system.addBasisSet(basisCoreA);
    system.addBasisSet(basisCoreB);

    system.setupOneParticleMatrix();
    system.setupTwoParticleMatrix();

    HFsolver solver(system);
    solver.runSolver();

    return 0;
}


