#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <armadillo>

#include<system/system.h>

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
    mat m_P;
    mat m_C;
    field<mat> m_Q;

    void setupTwoParticleMatrix();
    void setupDensityMatrix();
    void normalize();
};

#endif // HFSOLVER_H
