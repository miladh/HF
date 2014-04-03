#include "kineticintegralgd.h"

using namespace hf;

KineticIntegralGD::KineticIntegralGD(KineticIntegral* kinetic,
                                     OverlapIntegralGD *overlapGD):
    m_kinetic(kinetic),
    m_overlapGD(overlapGD),
    m_overlap(overlapGD->overlap()),
    m_primitiveA(m_overlap->primitiveA()),
    m_primitiveB(m_overlap->primitiveB())
{
}


double hf::KineticIntegralGD::evaluate(int cor, int iA, int iB)
{
    double b = m_primitiveB->exponent();

    double dS_iA_iBnn = m_overlapGD->evaluate(cor, iA, iB + 2);
    double dS_iA_iB   = m_overlapGD->evaluate(cor, iA, iB);
    double dS_iA_iBpp;

    if(iB - 2 >= 0) {
        dS_iA_iBpp= m_overlapGD->evaluate(cor, iA, iB - 2);
    } else {
        dS_iA_iBpp = 0;
    }
    return 4 * b * b * dS_iA_iBnn - 2*b * (2*iB + 1) * dS_iA_iB + iB * (iB - 1) * dS_iA_iBpp;
}

rowvec hf::KineticIntegralGD::evaluate()
{
    rowvec dT = zeros<rowvec>(3);

    double T_iA_iB = m_kinetic->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double T_jA_jB = m_kinetic->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double T_kA_kB = m_kinetic->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    double dT_iA_iB = evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double dT_jA_jB = evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double dT_kA_kB = evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    double S_iA_iB = m_overlap->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double S_jA_jB = m_overlap->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double S_kA_kB = m_overlap->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    double dS_iA_iB = m_overlapGD->evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower());
    double dS_jA_jB = m_overlapGD->evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower());
    double dS_kA_kB = m_overlapGD->evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower());

    dT(0) = (dT_iA_iB * S_jA_jB * S_kA_kB) + (dS_iA_iB * T_jA_jB * S_kA_kB) + (dS_iA_iB * S_jA_jB  * T_kA_kB);
    dT(1) = (T_iA_iB * dS_jA_jB * S_kA_kB) + (S_iA_iB * dT_jA_jB * S_kA_kB) + (S_iA_iB * dS_jA_jB  * T_kA_kB);
    dT(2) = (T_iA_iB * S_jA_jB * dS_kA_kB) + (S_iA_iB * T_jA_jB * dS_kA_kB) + (S_iA_iB * S_jA_jB  * dT_kA_kB);

    dT *= -0.5;
    return dT;

}
