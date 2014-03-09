#ifndef BOMD_H
#define BOMD_H

#include <iostream>
#include <armadillo>

#include "../system/system.h"
#include "../hfSolver/hfsolver.h"
#include "../hfSolver/rhf.h"
#include "../geometricalDerivative/geometricalderivative.h"

using namespace arma;
using namespace std;

namespace hf{

class BOMD
{
public:
    BOMD(System *system, HFsolver *solver, const int &rank, const int &nProcs);

    void runDynamics();
    void solveSingleStep();
    double getEnergy() const;

    const rowvec& getEnergyGradient() const;

private:
    int m_rank, m_nProcs;
    System* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    int m_nCores,  m_nOrbitals;
    int m_nSteps;

    double m_dtn;
    double m_dampingFactor;
    double m_energy, m_fockEnergy;


    rowvec m_energyGradient;
    mat pos, posNew, posOld;

    void IntegrateCoreForwardInTime(int core);
    void writeToFile(mat R, int currentTimeStep);
    void updateCorePositions();


};
}
#endif // BOMD_H


