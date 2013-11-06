#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <armadillo>

#include<integrator/integrator.h>
#include<basisSet/basisset.h>


using namespace arma;
using namespace std;

class System
{
public:
    System(int nElectrons, int maxAngularMomentum, rowvec coreCharges);

    void addBasisSet(BasisSet *basisSet);
    int getTotalNumOfBasisFunc();
    int getNumOfElectrons();

    double getNucleiPotential();
    rowvec getOneParticleIntegral(const int a, const int b);
    double getTwoParticleIntegral(const int p, const int q,
                                  const int r, const int s);

private:
    vector<BasisSet *> m_basisSet;
    Integrator integrator;
    vector<int > m_coreID;
    vector<int > m_cumSumContracted;

    rowvec m_coreCharges;
    int m_nElectrons;
};

#endif // SYSTEM_H
