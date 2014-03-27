#ifndef ELECTRONICSYSTEM_H
#define ELECTRONICSYSTEM_H

#include <iostream>
#include <armadillo>
#include <mpi.h>

#include "../integrator/integrator.h"
#include "../basisSet/basisset.h"
#include "../atom/atom.h"


using namespace arma;
using namespace std;

namespace hf
{

class ElectronicSystem
{
public:
    ElectronicSystem(const int& nSpinUpElectrons, const int& nSpinDownElectrons, const int& maxAngularMomentum);
    ElectronicSystem(const int& nElectrons, const int& maxAngularMomentum);

    void addBasisSet(BasisSet *basisSet);
    int getTotalNumOfBasisFunc();
    const int &nElectrons() const;
    const int &nSpinUpElectrons() const;
    const int &nSpinDownElectrons() const;

    int getNumOfCores();


    rowvec getNucleiPotential_derivative(int activeCore);
    mat getOneParticleDerivative(const int a, const int b, const int N);
    rowvec getTwoParticleIntegralDerivative(const int a, const int b, const int c, const int d,
                                            const int N);




    void addAtom(Atom *atom);
    double overlapIntegral(const int &p, const int &q);
    double oneParticleIntegral(const int &p, const int &q);
    double twoParticleIntegral(const int &p, const int &q,
                               const int &r, const int &s);
    double nuclearPotential();
    double gaussianProduct(const int &p, const int &q,
                           const double &x, const double &y, const double &z);


    void computePartialCharge(const mat &PS);

    rowvec nuclearPotentialGD(int activeCore);
private:
    Integrator integrator;
    vector<Atom *> m_atoms;
    vector<const ContractedGTO *> m_basisFunctions;

    int m_nElectrons;
    int m_nSpinUpElectrons;
    int m_nSpinDownElectrons;
    int m_nBasisFunctions;
};

}


#endif // ELECTRONICSYSTEM_H
