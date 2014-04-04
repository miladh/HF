#include "atom.h"

using namespace hf;

Atom::Atom(string basisFile, const rowvec& corePosition, const rowvec& coreVelocity)
{
    loadBasisFile(basisFile);
    setCorePosition(corePosition);
    setCoreVelocity(coreVelocity);
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

int Atom::nContractedGTOs() const
{
    return m_contractedGTOs.size();
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

    for(ContractedGTO &CGTO : m_contractedGTOs) {
            CGTO.setCenter(&m_corePosition);
        }

}
const rowvec &Atom::coreVelocity() const
{
    return m_coreVelocity;
}

void Atom::setCoreVelocity(const rowvec &coreVelocity)
{
    m_coreVelocity = coreVelocity;
}



void Atom::setCorePartialCharge(const double& corePartialCharge)
{
    m_corePartialCharge = corePartialCharge;
}








