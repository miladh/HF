#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <armadillo>

#include<primitiveGTO/primitiveGTO.h>
#include<contractedGTO/contractedGTO.h>
#include<integrator/integrator.h>
#include<basisSet/basisset.h>


using namespace arma;
using namespace std;

class System
{
public:
    System(int nOrbitals, int nNuclei, int maxAngularMomentum);

    void addBasisSet(BasisSet *basisSet);

    void setupOneParticleMatrix();
    void setupTwoParticleMatrix();

    mat getOneParticleMatrix() const;
    field<mat> getTwoParticleMatrix() const;
    mat getOverlapMatrix() const;







    void setupOneParticleMatrix2();
    void setupTwoParticleMatrix2();

private:
    mat m_h, m_S, m_R;
    field<mat> m_Q;
    vector<BasisSet *> m_basisSet;
    Integrator integrator;
    vector<int > m_coreID;
    vector<int > m_cumSumContracted;
};

#endif // SYSTEM_H
