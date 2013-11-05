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
    System(int nOrbitals, int maxAngularMomentum, rowvec coreCharges);

    void addBasisSet(BasisSet *basisSet);

    void setupOneParticleMatrix();
    void setupTwoParticleMatrix();

    mat getOneParticleMatrix() const;
    field<mat> getTwoParticleMatrix() const;
    mat getOverlapMatrix() const;



    void setupOneParticleMatrix2();
    void setupTwoParticleMatrix2();

    int getTotalNumOfBasisFunc();
    int getNumOfElectrons();

private:
    mat m_h, m_S;
    field<mat> m_Q;
    vector<BasisSet *> m_basisSet;
    Integrator integrator;
    vector<int > m_coreID;
    vector<int > m_cumSumContracted;

    rowvec m_coreCharges;
};

#endif // SYSTEM_H
