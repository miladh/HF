#ifndef ELECTRONICSYSTEM_H
#define ELECTRONICSYSTEM_H

#include <iostream>
#include <armadillo>
#include <mpi.h>

#include "../integrator/integrator.h"
#include "../atom/atom.h"


using namespace arma;
using namespace std;

namespace hf
{

class ElectronicSystem
{
public:
    ElectronicSystem(const int& maxAngularMomentum);

    void addAtom(Atom *atom);
    const int &nAtoms();
    const int &nElectrons() const;
    const int &nSpinUpElectrons() const;
    const int &nSpinDownElectrons() const;
    const int &nBasisFunctions();

    double nuclearPotential();
    double overlapIntegral(const int &p, const int &q);
    double oneParticleIntegral(const int &p, const int &q);
    double twoParticleIntegral(const int &p, const int &q,
                               const int &r, const int &s);



    rowvec nuclearPotentialGD(int activeCore);
    mat getOneParticleDerivative(const int a, const int b, const int N);
    rowvec getTwoParticleIntegralDerivative(const int a, const int b, const int c, const int d,
                                            const int N);


    vector<const ContractedGTO *> basisFunctions() const;
    vector<Atom *> m_atoms; //Should be private!!!

    vector<Atom *> atoms() const;

private:
    Integrator integrator;
    vector<const ContractedGTO *> m_basisFunctions;

    int m_nElectrons = 0;
    int m_nSpinUpElectrons = 0;
    int m_nSpinDownElectrons = 0;
    int m_nAtoms = 0;
    int m_nBasisFunctions = 0;

};

}


#endif // ELECTRONICSYSTEM_H
