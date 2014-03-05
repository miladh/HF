#ifndef HFSOLVER_H
#define HFSOLVER_H

#include <iostream>
#include <iomanip>
#include <armadillo>
#include <includes/defines.h>
#include <system/system.h>

using namespace arma;
using namespace std;

class HFsolver
{
public:
   HFsolver(System *system, const int &rank, const int &nProcs);
   void runSolver();

   field<mat> getQmatrix();
   mat gethmatrix();
   mat getSmatrix();

   mat getC() const;

   double getEnergy() const;
   double getFockEnergy() const;
   void setupTwoParticleMatrix();
   void setupOneParticleMatrix();

   mat getDensityMatrix() const;
   mat getF();


   int m_step;
private:
   int m_rank, m_nProcs;
   System *m_system;
   cube m_density;

   int m_nElectrons, m_nOrbitals;
   mat m_S, m_h, m_F, m_P, m_C;
   field<mat> m_Q;

   double m_energy, m_fockEnergy;

   void normalize();
   void solveSingle();
   void setupFockMatrix();
   void calculateEnergy();
   void calculateDensity();
   void densityOutput(const double &xMin, const double &xMax, const double &yMin, const double &yMax, const double &zMin, const double &zMax);
};

#endif // HFSOLVER_H
