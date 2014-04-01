#include "contractedGTO.h"

using namespace hf;
ContractedGTO::ContractedGTO()
{
}

void ContractedGTO::addPrimitive(PrimitiveGTO primitiveGTO)
{
    m_primitivesGTOs.push_back(primitiveGTO);
}

void ContractedGTO::setCenter(const rowvec& center)
{
    m_center = center;
}

const rowvec& ContractedGTO::center() const
{
    return m_center;
}

const vector<PrimitiveGTO> &ContractedGTO::primitivesGTOs() const
{
    return m_primitivesGTOs;
}

