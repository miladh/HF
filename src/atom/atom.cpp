#include "atom.h"

using namespace hf;

Atom::Atom(string basisFile, const rowvec& corePosition):
    m_corePosition(corePosition)
{
    loadBasisFile(basisFile);
}

void Atom::loadBasisFile(string fileName)
{
    TurbomoleParser parser(fileName);
    m_atomType   = parser.atomType();
    m_coreCharge = parser.atomCharge();
    m_coreMass   = parser.atomMass();
    m_angularMomentum = parser.angularMomentum();
    m_nElectrons = int(m_atomType);
    m_contractedGTOs = parser.contractedGTOs();

    cout << m_contractedGTOs.size() << endl;
    cout << m_angularMomentum << m_atomType  << m_nElectrons << m_coreCharge <<"   " << m_coreMass << endl;

}

const int& Atom::atomType() const
{
    return m_atomType;
}

const int& Atom::nElectrons() const
{
    return m_nElectrons;
}

const int& Atom::coreCharge() const
{
    return m_coreCharge;
}

const double& Atom::coreMass() const
{
    return m_coreMass;
}

const int& Atom::angularMomentum() const
{
    return m_angularMomentum;
}

const double& Atom::corePartialCharge() const
{
    return m_corePartialCharge;
}

const rowvec& Atom::corePosition() const
{
    return m_corePosition;
}

const vector<ContractedGTO>& Atom::contractedGTOs() const
{
    return m_contractedGTOs;
}

void Atom::setCorePosition(const rowvec &corePosition)
{
    m_corePosition = corePosition;
}

void Atom::setCorePartialCharge(const double& corePartialCharge)
{
    m_corePartialCharge = corePartialCharge;
}








