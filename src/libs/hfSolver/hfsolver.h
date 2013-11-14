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
    HFsolver(System *system);
    void runSolver();


    field<mat> getQmatrix();
    mat gethmatrix();
    mat getSmatrix();

private:
    System *m_system;

    int m_nElectrons, m_nOrbitals;
    mat m_S, m_h, m_F, m_P, m_C;
    field<mat> m_Q;
    double m_fockEnergy, m_energy, m_toler;

    void normalize();
    void solveSingle();
    void setupFockMatrix();
    void setupTwoParticleMatrix();
    void setupOneParticleMatrix();
};

#endif // HFSOLVER_H
