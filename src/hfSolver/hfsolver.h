#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "../defines.h"
#include "../system/system.h"

using namespace arma;
using namespace std;

namespace hf
{
class HFsolver
{
public:
   HFsolver(System *system, const int &rank, const int &nProcs);
   void runSolver();

   const field<mat> &getQmatrix();
   const mat& gethmatrix();
   const mat& getSmatrix();

   const double &getEnergy() const;
   void setupTwoParticleMatrix();
   void setupOneParticleMatrix();

protected:
   int m_rank, m_nProcs, m_step;
   System *m_system;
   cube m_density;

   int m_nElectrons, m_nBasisFunctions;
   int m_nSpinUpElectrons, m_nSpinDownElectrons;
   mat m_S, m_h;
   field<mat> m_Q;

   double m_energy;

   virtual void advance() = 0;
   virtual void normalize() = 0;
   virtual void solveSingle() = 0;
   virtual void updateFockMatrix() = 0;
   virtual void calculateEnergy()=0;
   virtual void calculateDensity()= 0;
   void densityOutput(const double &xMin, const double &xMax, const double &yMin, const double &yMax, const double &zMin, const double &zMax);
};
}
#endif // HFSOLVER_H
