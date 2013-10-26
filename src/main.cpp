#include <iostream>
#include <armadillo>
#include<src/system/system.h>
#include<src/integrator/integrator.h>
//#include<src/math/boys.h>

using namespace arma;
using namespace std;


int main()
{
    Integrator solver;
    mat pos = zeros(2,3);
    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    pos.row(0) = posA;
    pos.row(1) = posB;

    solver.setCorePositionA(posA);
    solver.setCorePositionB(posB);
    solver.setCorePositions(pos);
    solver.setMaxAngularMomentum(2);

    solver.addPrimitives(new PrimitiveGTO(0.2, 1.0));
    solver.addPrimitives(new PrimitiveGTO(0.3, 1.0));
    solver.setupE();
//    solver.setupR(posA);
//    Boys boys(2.3252,12);
//    rowvec Fn = boys.getBoysFunctions();
//    cout << scientific <<setprecision(10) << Fn[6] <<endl;
    return 0;
}

