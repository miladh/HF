#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>

#include<src/primitiveGTO/primitiveGTO.h>

using namespace arma;
using namespace std;

class Integrator
{
public:
    Integrator();

    int l = 2;
    int tmax;
    vector <PrimitiveGTO *> primitives;
    void setE();
    cube E[3];
    mat R = zeros<mat>(2,3);

    mat index = zeros<mat>(3,2);
};

#endif // INTEGRATOR_H
