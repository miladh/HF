#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <boost/mpi.hpp>
#include "../defines.h"
#include "../system/electronicsystem.h"


using namespace arma;
using namespace std;

namespace hf
{
class HFsolver
{
public:
   HFsolver(ElectronicSystem *system);
   void runSolver();

   const mat& overlapMatrix() const;
   virtual field<const mat *> fockMatrix() = 0;
   virtual field<const mat *> densityMatrix()  const= 0;
   virtual field<const mat *> expansionCoefficients()  const= 0;

   const double &energy() const;
   void setupTwoParticleMatrix();
   void setupOneParticleMatrix();



protected:
   ElectronicSystem * m_system;
   int m_rank;
   int m_nProcs;
   int m_iteration;

   int m_nElectrons;
   int m_nSpinUpElectrons;
   int m_nSpinDownElectrons;
   int m_nBasisFunctions;
   mat m_S;
   mat m_h;
   field<mat> m_Q;

   double m_energy;

   // MPI-----------------------
   boost::mpi::communicator m_world;
   boost::mpi::timer m_timer;
   ivec m_basisIndexToProcsMap;
   vector<int> m_myBasisIndices;
   //---------------------------

   virtual void advance() = 0;
   virtual void solveSingle() = 0;
   virtual void updateFockMatrix() = 0;
   virtual void calculateEnergy()=0;

   const mat &normalize(mat &C, const int &HOcoeff);
   double computeStdDeviation(const vec &fockEnergies, const vec &fockEnergiesOld);
};
}
#endif // HFSOLVER_H
