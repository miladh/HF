#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include "../includes/defines.h"
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

   field<mat> getQmatrix();
   mat gethmatrix();
   mat getSmatrix();

   double getEnergy() const;
   double getFockEnergy() const;
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

   double m_energy, m_fockEnergy;

   virtual void normalize() = 0;
   virtual void solveSingle() = 0;
   virtual void updateFockMatrix() = 0;
   virtual void calculateEnergy()=0;
   virtual void calculateDensity()= 0;
   void densityOutput(const double &xMin, const double &xMax, const double &yMin, const double &yMax, const double &zMin, const double &zMax);
};
}
#endif // HFSOLVER_H
