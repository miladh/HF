#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <mpi.h>
#include "../defines.h"
#include "../system/system.h"
#include "../system/electronicsystem.h"


using namespace arma;
using namespace std;

namespace hf
{
class HFsolver
{
public:
   HFsolver(ElectronicSystem *system, const int &rank, const int &nProcs);
   void runSolver();

   const field<mat> &getQmatrix();
   const mat& gethmatrix();
   const mat& getSmatrix();
   virtual field<mat> getFockMatrix() = 0;
   virtual field<mat> getDensityMatrix()  const= 0;

   const double &getEnergy() const;
   void setupTwoParticleMatrix();
   void setupOneParticleMatrix();



protected:
   int m_rank, m_nProcs, m_step, m_iteration;
   ElectronicSystem * m_system;
   cube m_density;

   int m_nElectrons, m_nSpinUpElectrons, m_nSpinDownElectrons, m_nBasisFunctions;
   mat m_S, m_h;
   field<mat> m_Q;

   double m_energy;

   // MPI-----------------------
   ivec m_basisIndexToProcsMap;
   vector<int> m_myBasisIndices;
   //---------------------------

   virtual void advance() = 0;
   virtual void solveSingle() = 0;
   virtual void updateFockMatrix() = 0;
   virtual void calculateEnergy()=0;
   virtual void calculateDensity()= 0;


   void densityOutput(const double &xMin, const double &xMax, const double &yMin, const double &yMax, const double &zMin, const double &zMax);
   const mat &normalize(mat &C, const int &HOcoeff);
   double computeStdDeviation(const vec &fockEnergies, const vec &fockEnergiesOld);
};
}
#endif // HFSOLVER_H
