#include <iostream>
#include <armadillo>

#include<src/system/system.h>
#include<src/hfSolver/hfsolver.h>
#include<src/contractedGTO/contractedGTO.h>
#include<src/basisSet/quadzeta.h>


using namespace arma;
using namespace std;


int main()
{

    int nBasisFunc = 4;
    int nNuclei    = 2;
    int nSteps     = 20;
    int maxAngularMomentum = 0;

    int nOrbitals = nBasisFunc * nNuclei;
    System system(nOrbitals, nNuclei, maxAngularMomentum);


    QuadZeta *basisCoreA  = new QuadZeta;
    QuadZeta *basisCoreB  = new QuadZeta;

    rowvec A = {-0.5 , 0, 0};
    rowvec B = {0.5 , 0, 0};
    basisCoreA->setCorePosition(A);
    basisCoreB->setCorePosition(B);

    system.addBasisSet(basisCoreA);
    system.addBasisSet(basisCoreB);

    system.setupOneParticleMatrix();
    system.setupTwoParticleMatrix();

    HFsolver solver(system,nOrbitals,nSteps);
    solver.runSolver();

    return 0;
}

