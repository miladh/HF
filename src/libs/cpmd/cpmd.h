#ifndef CPMD_H
#define CPMD_H

#include <iostream>
#include <iomanip>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/splitValence/h_321g.h>
#include<basisSet/splitValence/h_431g.h>
#include<basisSet/splitValence/li_321g.h>
#include<basisSet/splitValence/o_321g.h>
#include<basisSet/splitValence/o_431g.h>

using namespace arma;
using namespace std;

class CPMD
{
public:
    CPMD(System *system);

    void runDynamics();

private:
    System* m_system;
    HFsolver *m_solver;

    int m_nCores, m_nElectrons, m_nOrbitals;
    int m_nSteps, m_eSteps;

    double m_dte, m_dtn;
    double m_gammaE, m_gammaN;
    double m_massE, m_energy;


    rowvec m_lambda, m_energyGradient;
    mat m_C, m_Cp, m_Cm;
    mat m_S, m_h, m_F, m_P;
    field<mat> m_Q;

    mat pos, posNew, posOld;

    field<field<rowvec>> m_dQ;
    field<rowvec> m_dS, m_dh;


    void setupFockMatrix();
    void setupDerivativeMatrices(const int core);
    void IntegrateWavefunctionForwardInTime(int orb);
    void IntegrateCoreForwardInTime(int core);
    double calculateEnergy();
    rowvec calculateEnergy_derivative(int core);
    mat normalize(mat C, mat S);
    void writeToFile(const mat R, int n);
    void updateCorePositions();
};

#endif // CPMD_H
