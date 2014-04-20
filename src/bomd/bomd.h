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
    void computeForces();

    const mat &energyGradient() const;
    double potentialEnergy() const;

private:
    ElectronicSystem* m_system;
    HFsolver *m_solver;
    GeometricalDerivative* m_GD;
    vector<Atom *> m_atoms;

    int m_nAtoms;
    int m_nSteps;
    int m_rank;

    double m_dt;
    double m_frictionConstant;

    mat m_energyGradient;

    vec m_time;
    vec m_totalEnergy;
    vec m_kineticEnergy;
    vec m_potentialEnergy;


    void solveSingleStep();
    void initialStep();
    void halfKick();
    void updateCores();

    void writeLammpsFile(int currentTimeStep);
    void writeSystemProperties();


    void systemProperties(int currentTimeStep);
};
}
#endif // BOMD_H


