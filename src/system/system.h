#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <armadillo>
#include <mpi.h>

#include "../integrator/integrator.h"
#include "../basisSet/basisset.h"


using namespace arma;
using namespace std;

namespace hf
{

class System
{
public:
    System(const int& nSpinUpElectrons, const int& nSpinDownElectrons, const int& maxAngularMomentum);
    System(const int& nElectrons, const int& maxAngularMomentum);

    vector<BasisSet *> m_basisSet;

    void addBasisSet(BasisSet *basisSet);
    int getTotalNumOfBasisFunc();
    const int &getNumOfElectrons();
    const int &getNumOfSpinUpElectrons();
    const int &getNumOfSpinDownElectrons();

    int getNumOfCores();

    double getNucleiPotential();
    rowvec getOneParticleIntegral(const int a, const int b);
    double getTwoParticleIntegral(const int p, const int q,
                                  const int r, const int s);

    rowvec getNucleiPotential_derivative(int activeCore);
    mat getOneParticleDerivative(const int a, const int b, const int N);
    rowvec getTwoParticleIntegralDerivative(const int a, const int b, const int c, const int d,
                                            const int N);


    mat getOneParticleDerivativeOfOrbitals(const int a, const int b, const int N);
    rowvec getTwoParticleIntegralDerivativeOfOrbitals(const int a, const int b, const int c, const int d, const int N);
    double gaussianProduct(const int a, const int b, const double &x, const double &y, const double &z);

    void computePartialCharge(const mat &PS);
private:
    Integrator integrator;
    vector<int > m_coreID;
    vector<int > m_cumSumContracted;

    int m_nElectrons, m_nSpinUpElectrons,m_nSpinDownElectrons;
};

}
#endif // SYSTEM_H
