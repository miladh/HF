#ifndef GEOMETRICALDERIVATIVE_H
#define GEOMETRICALDERIVATIVE_H


#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include "../hfSolver/hfsolver.h"
#include "../hfSolver/rhf.h"

using namespace arma;
using namespace std;

namespace hf {

class GeometricalDerivative
{
public:
    GeometricalDerivative(ElectronicSystem *system, HFsolver *solver);
    const mat &energyGradient();

private:
    ElectronicSystem *m_system;
    HFsolver *m_solver;
    mat m_gradE, m_totGradE;

    int m_nBasisFunctions;
    int m_rank, m_nProcs;

    void calculateEnergyGradient();

    // MPI-----------------------
    boost::mpi::communicator m_world;
    boost::mpi::timer m_timer;
    imat m_pqIndicesToProcsMap;
    vector<pair<int,int> > m_myPQIndices;
    //---------------------------
};

}
#endif // GEOMETRICALDERIVATIVE_H
