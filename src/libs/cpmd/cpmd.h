#ifndef CPMD_H
#define CPMD_H

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

class cpmd
{
public:
    cpmd();

    void runDynamics();

private:
    int m_nElectrons, m_nOrbitals;
    int m_nSteps, m_eSteps;

    double m_dte, m_dtn;
    double m_gammaE, m_gammaN;
    double m_massE, m_massN;
    double m_lambda;

    rowvec m_energyGradient;
    rowvec m_coreCharges;
    vec m_C, m_Cp, m_Cm;
    mat m_S, m_h, m_F;

    mat pos, posNew, posOld;

    field<mat> m_Q, m_dQ;
    field<rowvec> m_dS, m_dh;

    BasisSet *m_basisCoreA;
    BasisSet *m_basisCoreB;
    System* m_system;
    HFsolver *m_solver;



    void systemConfiguration();
    void setupFockMatrix();
    void setupDerivativeMatrices(const int core);
    void IntegrateWavefunctionForwardInTime();
    void IntegrateCoreForwardInTime(int core);
    double calculateEnergy(field<mat> Q, const mat h, const vec C);
    rowvec calculateEnergy_derivative(const field<mat> dQ, const field<rowvec> dh, const vec C, int core);
    vec normalize(vec C, mat S);
    void writeToFile(const mat R, int n);
};

#endif // CPMD_H
