#include "contractedGTO.h"

ContractedGTO::ContractedGTO()
{
}

vector<PrimitiveGTO*> ContractedGTO::primitives() const
{
    return m_primitives;
}

void ContractedGTO::addPrimitiveGTO(PrimitiveGTO *primitive)
{
    m_primitives.push_back(primitive);
}
