#include "overlapintegralgd.h"

using namespace hf;


OverlapIntegralGD::OverlapIntegralGD(OverlapIntegral *overlap,
                                     const field<cube> *dEab_dQab):
    m_overlap(overlap),
    m_dEab_dQab(dEab_dQab),
    m_primitiveA(overlap->primitiveA()),
    m_primitiveB(overlap->primitiveB())

{
}

double OverlapIntegralGD::evaluate(int cor, int iA, int iB)
{
     return m_dEab_dQab->at(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA->exponent() + m_primitiveB->exponent()));
}

rowvec OverlapIntegralGD::evaluate()
{
    rowvec dSab_dQab = zeros<rowvec>(3);

    dSab_dQab(0) =   evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower())
                    * m_overlap->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower())
                    * m_overlap->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    dSab_dQab(1) =   evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower())
                    * m_overlap->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower())
                    * m_overlap->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    dSab_dQab(2) =   evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower())
                    * m_overlap->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower())
                    * m_overlap->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());

    return dSab_dQab * m_primitiveA->weight() * m_primitiveB->weight();
}

OverlapIntegral *OverlapIntegralGD::overlap() const
{
    return m_overlap;
}


