#ifndef ANALYSER_H
#define ANALYSER_H

#include <armadillo>
#include <iostream>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include "../hfSolver/hfsolver.h"
#include "../system/electronicsystem.h"
#include "../outputManager/outputmanager.h"


using namespace arma;
using namespace std;

namespace hf {

class Analyser
{
public:
    Analyser(ElectronicSystem *system, HFsolver *solver);

    void saveEnergies();
    void atomicPartialCharge();
    void calculateChargeDensity();
    void calculateElectrostaticPotential();
    void dipoleMoment();
    void saveResults();



private:
    ElectronicSystem* m_system;
    HFsolver* m_solver;
    Integrator* m_integrator;
    OutputManager *m_outputManager;
    vector<const ContractedGTO *> m_basisFunctions;

    int m_nBasisFunctions;
    int m_rank;
    int m_nProcs;


    double electronicPotential(const int &p, const int &q, const rowvec &P);
    void writeDensityToFile(const cube &density, const double &xMin, const double &xMax,
                       const double &yMin, const double &yMax,
                       const double &zMin, const double &zMax);

    double gaussianProduct(const int &p, const int &q,
                           const double &x, const double &y, const double &z);
};


}
#endif // ANALYSER_H
