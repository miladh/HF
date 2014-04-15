#include "dipoleintegral.h"

using namespace hf;

DipoleIntegral::DipoleIntegral(const field<cube>* Eab,
                               const PrimitiveGTO* primitiveA,
                               const PrimitiveGTO* primitiveB):
  m_Eab(Eab),
  m_primitiveA(primitiveA),
  m_primitiveB(primitiveB)
{
}


double DipoleIntegral::evaluate(int cor, int iA, int iB, const double& P)
{

    return  (m_Eab->at(cor)(iA,iB,1) + P * m_Eab->at(cor)(iA,iB,0))
            * sqrt(M_PI / (m_primitiveA->exponent() + m_primitiveB->exponent()));
}

double DipoleIntegral::evaluate(int cor, int iA, int iB)
{
    return m_Eab->at(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA->exponent() + m_primitiveB->exponent()));
}

rowvec DipoleIntegral::evaluate()
{
    const rowvec &A = m_primitiveA->center();
    const rowvec &B = m_primitiveB->center();

    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    double p = a + b;
    rowvec P = (a*A + b*B)/p;
    rowvec D = {0,0,0};

    double S_iA_iB = evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double S_jA_jB = evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double S_kA_kB = evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    D(0) = evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower(), P(0))
         * S_jA_jB * S_kA_kB;


    D(1) = evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower(), P(1))
            * S_iA_iB * S_kA_kB;

    D(2) = evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower(), P(2))
            * S_jA_jB * S_iA_iB;

    D *= m_primitiveA->weight() * m_primitiveB->weight();

    return D;
}


const PrimitiveGTO *DipoleIntegral::primitiveA() const
{
    return m_primitiveA;
}
const PrimitiveGTO *DipoleIntegral::primitiveB() const
{
    return m_primitiveB;
}
