#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include<includes/defines.h>
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

    mat getC() const;

    double getEnergy() const;

private:
    System *m_system;

    int m_nElectrons, m_nOrbitals;
    mat m_S, m_h, m_F, m_P, m_C;
    field<mat> m_Q;
    double m_energy, m_fockEnergy;

    void normalize();
    void solveSingle();
    void setupFockMatrix();
    void setupTwoParticleMatrix();
    void setupOneParticleMatrix();
    void calculateEnergy();
};

#endif // HFSOLVER_H
