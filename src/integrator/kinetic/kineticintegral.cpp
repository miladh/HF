#include "kineticintegral.h"

using namespace hf;
KineticIntegral::KineticIntegral(OverlapIntegral *overlap,
                                 const PrimitiveGTO* primitiveA,
                                 const PrimitiveGTO* primitiveB):

    m_overlap(overlap),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB)
{
}

double KineticIntegral::evaluate(int cor, int iA, int iB)
{
    double b = m_primitiveB->exponent();

    double S_iA_iBnn = m_overlap->evaluate(cor, iA, iB + 2);
    double S_iA_iB = m_overlap->evaluate(cor, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= m_overlap->evaluate(cor, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB - 1) * S_iA_iBpp;

}

double KineticIntegral::evaluate()
{
    double T_iA_iB = evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double T_jA_jB = evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double T_kA_kB = evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    double S_iA_iB = m_overlap->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double S_jA_jB = m_overlap->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double S_kA_kB = m_overlap->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5 * m_primitiveA->weight() * m_primitiveB->weight();
    return result;

}
