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

    for(PrimitiveGTO &PGTO : m_primitivesGTOs){
        PGTO.setCenter(m_center);
    }
}

const rowvec& ContractedGTO::center() const
{
    return m_center;
}

const vector<PrimitiveGTO> &ContractedGTO::primitiveGTOs() const
{
    return m_primitivesGTOs;
}

