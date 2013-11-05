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
    HFsolver(System system);
    void runSolver();

private:
    System m_system;
    int m_nElectrons,m_nOrbitals;
    double m_fockEnergy, m_energy, m_toler;
    mat m_F, m_S, m_h, m_P, m_C;
    field<mat> m_Q;

    void setupFockMatrix();
    void normalize();
    void solveSingle();
};

#endif // HFSOLVER_H
