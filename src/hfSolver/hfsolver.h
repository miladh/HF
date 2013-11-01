#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <armadillo>

#include<src/system/system.h>

using namespace arma;
using namespace std;

class HFsolver
{
public:
    HFsolver(System system, int nOrbitals, int nSteps);
    void runSolver();

private:
    System m_system;
    int m_nSteps;

    mat m_F, m_S, m_G, m_h;
    vec m_C;
    field<mat> m_Q;

    void setupTwoParticleMatrix();
    void normalize();
};

#endif // HFSOLVER_H
