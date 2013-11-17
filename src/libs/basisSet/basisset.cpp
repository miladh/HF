#include "basisset.h"

BasisSet::BasisSet()
{
}
rowvec BasisSet::corePosition() const
{
    return m_corePosition;
}

void BasisSet::setCorePosition(const rowvec &corePosition)
{
    m_corePosition = corePosition;
}

const ContractedGTO &BasisSet::getContracted(const int c) const
{
    return m_contractedGTOs.at(c);
}

int BasisSet::getNumContracted() const
{
    return m_contractedGTOs.size();
}
int BasisSet::coreCharge() const
{
    return m_coreCharge;
}

void BasisSet::setCoreCharge(const int &coreCharge)
{
    m_coreCharge = coreCharge;
}
int BasisSet::coreMass() const
{
    return m_coreMass;
}

void BasisSet::setCoreMass(const int &coreMass)
{
    m_coreMass = PROTONMASS * coreMass;
}





