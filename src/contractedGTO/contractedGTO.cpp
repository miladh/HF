#include "contractedGTO.h"

using namespace hf;
ContractedGTO::ContractedGTO()
{
}

void ContractedGTO::addPrimitive(PrimitiveGTO primitiveGTO)
{
    m_primitivesGTOs.push_back(primitiveGTO);
}

const PrimitiveGTO& ContractedGTO::getPrimitive(const int p) const
{
    return m_primitivesGTOs.at(p);
}

int ContractedGTO::getNumPrimitives() const
{
    return m_primitivesGTOs.size();
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

void ContractedGTO::setPrimitivesGTOs(const vector<PrimitiveGTO> &primitivesGTOs)
{
    m_primitivesGTOs = primitivesGTOs;
}
