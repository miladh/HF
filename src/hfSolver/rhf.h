#ifndef RHF_H
#define RHF_H

#include <armadillo>
#include <iostream>
#include "../hfSolver/hfsolver.h"

using namespace arma;
using namespace std;

namespace hf {

class RHF : public HFsolver
{
public:
    RHF(System *system, const int &rank, const int &nProcs);

    mat getFockMatrix();
    mat getDensityMatrix() const;
    mat getExpansionCoeff() const;

private:
    mat m_F, m_C, m_P;
    int m_nElectrons;
    void solveSingle();
    void setupFockMatrix();
    void calculateEnergy();
    void updateFockMatrix();
    void calculateDensity();
    void normalize();
};

}

#endif // RHF_H
