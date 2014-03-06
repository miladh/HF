#ifndef BOMD_H
#define BOMD_H

#include <iostream>
#include <armadillo>

#include <system/system.h>
#include <hfSolver/hfsolver.h>
#include <geometricalDerivative/geometricalderivative.h>

using namespace arma;
using namespace std;
using namespace hf;

class BOMD
{
public:
    BOMD(System *system, const int &rank, const int &nProcs);

    void runDynamics();
    void solveSingleStep();
    rowvec getEnergyGradient() const;
    double getEnergy() const;

private:
    int m_rank, m_nProcs;
    System* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    int m_nCores, m_nElectrons, m_nOrbitals;
    int m_nSteps;

    double m_dtn;
    double m_dampingFactor;
    double m_energy, m_fockEnergy;


    rowvec m_energyGradient, m_pulayForce;
    mat m_S, m_h, m_P;

    mat pos, posNew, posOld;

    field<mat> m_Q;
    field<field<rowvec>>m_dQ;
    field<rowvec> m_dS, m_dh;
    field<rowvec> m_pulayS, m_pulayh;
    field<field<rowvec>>m_pulayQ;



    void IntegrateWavefunctionForwardInTime(int orb);
    void IntegrateCoreForwardInTime(int core);

    void writeToFile(mat R, int currentTimeStep);
    void updateCorePositions();
    void setupPulayMatrices(const int core);
    rowvec calculatePulayForce(int core);

};

#endif // BOMD_H


