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
    BOMD(System *system);
    void runDynamics();


private:
    System* m_system;
    HFsolver *m_solver;
    int m_nCores, m_nElectrons, m_nOrbitals;
    int m_nSteps;

    double m_dtn;
    double m_dampingFactor;
    double m_energy;


    rowvec m_energyGradient;
    mat m_C, m_Cp, m_Cm;
    mat m_S, m_h, m_P;

    mat pos, posNew, posOld;

    field<mat> m_Q;
    field<field<rowvec>>m_dQ;
    field<rowvec> m_dS, m_dh;


    void setupDerivativeMatrices(const int core);

    void IntegrateWavefunctionForwardInTime(int orb);
    void IntegrateCoreForwardInTime(int core);

    rowvec calculateEnergyGradient(int core);
    void writeToFile(const mat R, int n);
    void updateCorePositions();
};

#endif // BOMD_H


