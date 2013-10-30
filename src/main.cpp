#include <iostream>
#include <armadillo>
#include<src/system/system.h>
#include<src/integrator/integrator.h>


using namespace arma;
using namespace std;


int main()
{
    Integrator solver;
    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec posC = {2.3,0.9,3.2};
    rowvec posD = {-2.0,1.9,2.2};

    solver.setCorePositionA(posA);
    solver.setCorePositionB(posB);
    solver.setCorePositionC(posC);
    solver.setCorePositionD(posD);
    solver.setMaxAngularMomentum(3);

    solver.addPrimitives(new PrimitiveGTO(0.2, 1.0));
    solver.addPrimitives(new PrimitiveGTO(0.3, 1.0));
    solver.addPrimitives(new PrimitiveGTO(0.4, 1.0));
    solver.addPrimitives(new PrimitiveGTO(0.5, 1.0));
    solver.setupE();


//    double Ven = solver.nuclearAttractionIntegral(2,1,0,1,1,0);
//    cout << scientific << setprecision(10) << Ven << endl;

    double Vee = solver.electronRepulsionIntegral(1,2,0,1,0,0,0,0,0,0,0,3);
    cout << Vee <<endl;

    return 0;
}

