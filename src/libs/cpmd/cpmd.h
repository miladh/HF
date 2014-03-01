#ifndef CPMD_H
#define CPMD_H

#include <iostream>
#include <iomanip>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>

using namespace arma;
using namespace std;

class CPMD
{
public:
    CPMD(System *system, const int &rank, const int &nProcs);

    void runDynamics();

private:
    int m_rank, m_nProcs;
    System* m_system;
    HFsolver *m_solver;

    int m_nCores, m_nElectrons, m_nOrbitals;
    int m_nSteps, m_eSteps;

    double m_dte, m_dtn;
    double m_gammaE, m_gammaN;
    double m_massE, m_energy;


    mat m_C, m_Cp, m_Cm;
    mat m_S, m_h, m_F, m_P;
    mat m_lambda;
    mat pos, posOld, posNew;

    rowvec m_energyGradient;

    field<mat> m_Q;
    field<field<rowvec>> m_dQ;
    field<rowvec> m_dS, m_dh;


    void setupFockMatrix();
    void setupDerivativeMatrices(const int core);
    void IntegrateWavefunctionForwardInTime(int orb);
    void IntegrateCoreForwardInTime(int core);
    void setupLambdaMatrix(int orb);
    void updateCorePositions();
    void writeToFile(const mat R, int n);
    double calculateEnergy();
    rowvec calculateEnergy_derivative(int core);
    mat normalize(mat C, mat S);
};

#endif // CPMD_H
