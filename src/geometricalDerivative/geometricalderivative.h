#ifndef GEOMETRICALDERIVATIVE_H
#define GEOMETRICALDERIVATIVE_H


#include <iostream>
#include <armadillo>
#include "../system/system.h"
#include "../hfSolver/hfsolver.h"
#include "../hfSolver/rhf.h"

using namespace arma;
using namespace std;

namespace hf {

class GeometricalDerivative
{
public:
    GeometricalDerivative(System *system, HFsolver *solver);
    const rowvec &energyGradient(const int core);

private:
    System *m_system;
    HFsolver *m_solver;
    rowvec3 m_gradE;
    field<field<rowvec3>>m_dQ;
    field<rowvec3> m_dS, m_dh;

    int m_nBasisFunctions;
    int m_differentiationCore;

    void setupDerivativeMatrices();
    void calculateEnergyGradient();
};

}
#endif // GEOMETRICALDERIVATIVE_H
