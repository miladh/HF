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

    const mat& getFockMatrix();
    const mat& getDensityMatrix() const;
    const mat& getExpansionCoeff() const;

private:
    mat m_F, m_C, m_P;
    int m_nElectrons;
    double m_fockEnergy;

protected:
    void advance();
    void solveSingle();
    void calculateEnergy();
    void updateFockMatrix();
    void calculateDensity();
    void normalize();
};

}

#endif // RHF_H
