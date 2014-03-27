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
    ElectronicSystem(const int& maxAngularMomentum);


    vector<Atom *> m_atoms; //Should be private!!!

    void addBasisSet(BasisSet *basisSet);
    int nBasisFunctions();
    int nAtoms();
    const int &nElectrons() const;
    const int &nSpinUpElectrons() const;
    const int &nSpinDownElectrons() const;

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
    vector<const ContractedGTO *> m_basisFunctions;

    int m_nElectrons = 0;
    int m_nSpinUpElectrons = 0;
    int m_nSpinDownElectrons = 0;

};

}


#endif // ELECTRONICSYSTEM_H
