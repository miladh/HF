#include "overlapintegral.h"

using namespace hf;

OverlapIntegral::OverlapIntegral(const field<cube>* Eab,
                                 const PrimitiveGTO* primitiveA,
                                 const PrimitiveGTO* primitiveB):
    m_Eab(Eab),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB)
{
}


double OverlapIntegral::evaluate(int cor, int iA, int iB)
{
    return m_Eab->at(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA->exponent() + m_primitiveB->exponent()));
}

double OverlapIntegral::evaluate()
{

    return    evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower())
            * evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower())
            * evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower())
            * m_primitiveA->weight() * m_primitiveB->weight();
}
const PrimitiveGTO *OverlapIntegral::primitiveA() const
{
    return m_primitiveA;
}
const PrimitiveGTO *OverlapIntegral::primitiveB() const
{
    return m_primitiveB;
}



