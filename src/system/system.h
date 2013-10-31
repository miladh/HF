#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <armadillo>

#include<src/primitiveGTO/primitiveGTO.h>
#include<src/contractedGTO/contractedGTO.h>
#include<src/integrator/integrator.h>

using namespace arma;
using namespace std;

class System
{
public:
    System(int nOrbitals, int nNuclei, int maxAngularMomentum);
    void addPrimitives(PrimitiveGTO *primitive);


    void setupOneParticleMatrix();
    void setupTwoParticleMatrix();

    mat getOneParticleMatrix() const;
    double ****getTwoParticleMatrix() const;
    mat getOverlapMatrix() const;
private:
    mat m_h, m_S, m_R;
    double ****m_Q;
    vector<PrimitiveGTO *> m_primitives;
    Integrator integrator;
};

#endif // SYSTEM_H
