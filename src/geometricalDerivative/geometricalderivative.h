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

    void calculateEnergyGradient();

    int m_rank, m_nProcs;
    // MPI-----------------------
    boost::mpi::communicator m_world;
    boost::mpi::timer m_timer;
    ivec m_basisIndexToProcsMap;
    vector<int> m_myBasisIndices;
    //---------------------------
};

}
#endif // GEOMETRICALDERIVATIVE_H
