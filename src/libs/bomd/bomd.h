#ifndef BOMD_H
#define BOMD_H

#include <iostream>
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

class BOMD
{
public:
    BOMD();
    void runDynamics();

private:
    int m_nElectrons, m_nOrbitals;
    int m_nSteps;

    double m_dtn;
    double m_gammaN;
    double m_massN;
    double m_lambda;

    rowvec m_energyGradient;
    rowvec m_coreCharges;
    mat m_C, m_Cp, m_Cm;
    mat m_S, m_h, m_P;

    mat pos, posNew, posOld;

    field<mat> m_Q, m_dQ;
    field<rowvec> m_dS, m_dh;

    BasisSet *m_basisCoreA;
    BasisSet *m_basisCoreB;
    BasisSet *m_basisCoreC;
    System* m_system;
    HFsolver *m_solver;



    void systemConfiguration();
    void setupDerivativeMatrices(const int core);

    void IntegrateWavefunctionForwardInTime(int orb);
    void IntegrateCoreForwardInTime(int core);

    rowvec calculateEnergy_derivative(int core);
    void writeToFile(const mat R, int n);
};

#endif // BOMD_H


