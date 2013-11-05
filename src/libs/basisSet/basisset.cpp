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

vector<ContractedGTO *> BasisSet::contractedGTOs() const
{
    return m_contractedGTOs;
}


int BasisSet::getNumContracted() const
{
    return m_contractedGTOs.size();
}



