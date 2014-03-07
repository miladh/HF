#ifndef GEOMETRICALDERIVATIVE_H
#define GEOMETRICALDERIVATIVE_H


#include <iostream>
#include <armadillo>
#include <system/system.h>
#include <hfSolver/hfsolver.h>
#include <hfSolver/rhf.h>

using namespace arma;
using namespace std;

namespace hf {

class GeometricalDerivative
{
public:
    GeometricalDerivative(System *system, RHF *solver);
    rowvec energyGradient(const int core);

private:
    System *m_system;
    RHF *m_solver;
    field<field<rowvec>>m_dQ;
    field<rowvec> m_dS, m_dh;

    int m_nBasisFunctions;
    int m_differentiationCore;

    void setupDerivativeMatrices();
    rowvec calculateEnergyGradient();
};

}
#endif // GEOMETRICALDERIVATIVE_H
