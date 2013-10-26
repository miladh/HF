#include <iostream>
#include <armadillo>
#include<src/system/system.h>
#include<src/integrator/integrator.h>
//#include<src/math/boys.h>

using namespace arma;
using namespace std;


int main()
{
//    Integrator solver;
//    mat pos = zeros(2,3);
//    rowvec posA = {1.2,2.3,3.4};
//    rowvec posB = {-1.3,1.4,-2.4};
//    pos.row(0) = posA;
//    pos.row(1) = posB;

//    solver.setCorePositionA(posA);
//    solver.setCorePositionB(posB);
//    solver.setCorePositions(pos);
//    solver.setMaxAngularMomentum(2);

//    solver.addPrimitives(new PrimitiveGTO(0.2, 1.0));
//    solver.addPrimitives(new PrimitiveGTO(0.3, 1.0));
//    solver.setupE();

//    Boys boysF0_small(0);
//    mat F0; rowvec xvec;
//    xvec = linspace<rowvec>(0.01,50,20);
//    F0 = zeros(20,1);
//    for(uint i = 0; i < 20; i++){
//        boysF0_small.evaluateBoysFunctions(xvec[i]);
//        F0.row(i) =  boysF0_small.getBoysFunctions();
//    }

    return 0;
}

