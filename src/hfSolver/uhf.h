#ifndef UHF_H
#define UHF_H

#include <armadillo>
#include <iostream>
#include "../hfSolver/hfsolver.h"

using namespace arma;
using namespace std;

namespace hf {


class UHF : public HFsolver
{
public:
    UHF(System *system, const int &rank, const int &nProcs);

private:
    mat m_Fu, m_Cu, m_Pu;
    mat m_Fd, m_Cd, m_Pd;

    void solveSingle();
    void setupFockMatrix();
    void calculateEnergy();
    void updateFockMatrix();
    void calculateDensity();
    void normalize();
};

}
#endif // UHF_H
