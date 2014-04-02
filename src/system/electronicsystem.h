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
    ElectronicSystem();

    void addAtoms(vector<Atom *> atoms);
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


    vector<Atom *> atoms() const;
    vector<const ContractedGTO *> basisFunctions() const;

    mat overlapIntegralGD(const int &q, const int &p);
private:
    Integrator integrator;
    vector<Atom *> m_atoms;
    vector<const ContractedGTO *> m_basisFunctions;
    vector<int> m_basisFucntionIndexToAtomID;
    int m_nElectrons = 0;
    int m_nSpinUpElectrons = 0;
    int m_nSpinDownElectrons = 0;
    int m_nAtoms = 0;
    int m_nBasisFunctions = 0;

};

}


#endif // ELECTRONICSYSTEM_H
