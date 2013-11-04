#include <iostream>
#include <armadillo>

#include<src/system/system.h>
#include<src/hfSolver/hfsolver.h>


using namespace arma;
using namespace std;


int main()
{

    int nBasisFunc = 4;
    int nNuclei    = 2;
    int nSteps     = 20;
    int maxAngularMomentum = 0;

    int nOrbitals = nBasisFunc * nNuclei;
    vec exponent = {13.00773, 1.962079, 0.444529, 0.1219492};
    System system(nOrbitals, nNuclei, maxAngularMomentum);

    system.addPrimitives(new PrimitiveGTO(exponent[0],1.0));
    system.addPrimitives(new PrimitiveGTO(exponent[1],1.0));
    system.addPrimitives(new PrimitiveGTO(exponent[2],1.0));
    system.addPrimitives(new PrimitiveGTO(exponent[3],1.0));

    system.setupOneParticleMatrix();
    system.setupTwoParticleMatrix();

    HFsolver solver(system,nOrbitals,nSteps);
    solver.runSolver();

    return 0;
}

