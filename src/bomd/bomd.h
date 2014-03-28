#ifndef BOMD_H
#define BOMD_H

#include <iostream>
#include <armadillo>

#include "../hfSolver/hfsolver.h"
#include "../hfSolver/rhf.h"
#include "../geometricalDerivative/geometricalderivative.h"

using namespace arma;
using namespace std;

namespace hf{

class BOMD
{
public:
    BOMD(ElectronicSystem *system, HFsolver *solver);

    void runDynamics();
    void solveSingleStep();
    double getEnergy() const;

    const rowvec& getEnergyGradient() const;

private:
    int m_rank;
    ElectronicSystem* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    vector<Atom *> m_atoms;

    int m_nAtoms;
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


