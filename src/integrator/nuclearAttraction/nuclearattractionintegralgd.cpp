#include "nuclearattractionintegralgd.h"


using namespace hf;

NuclearAttractionIntegralGD::NuclearAttractionIntegralGD(const int highestOrder,
                                                         const field<cube> *Eab,
                                                         const field<cube> *dEab_dQab,
                                                         const PrimitiveGTO *primitiveA,
                                                         const PrimitiveGTO *primitiveB,
                                                         const rowvec *sourceCharge):
    m_R(new HermiteIntegrals(highestOrder)),
    m_Eab(Eab),
    m_dEab_dQab(dEab_dQab),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB),
    m_sourceCharge(sourceCharge)
{
}

void NuclearAttractionIntegralGD::updateHermiteIntegrals()
{
    const rowvec &A = m_primitiveA->center();
    const rowvec &B = m_primitiveB->center();
    const rowvec &C = (*m_sourceCharge);

    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    m_p = a + b;
    rowvec PC = (a*A + b*B)/m_p - C;

    m_tMax = m_primitiveA->xPower() + m_primitiveB->xPower() + 1;
    m_uMax = m_primitiveA->yPower() + m_primitiveB->yPower() + 1;
    m_vMax = m_primitiveA->zPower() + m_primitiveB->zPower() + 1;

    m_R->updateR(PC, m_p, m_tMax, m_uMax, m_vMax);
}



rowvec NuclearAttractionIntegralGD::QDerivative()
{
    rowvec dVdQ = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();


    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){
                dVdQ(0) += m_dEab_dQab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u)  * m_Eab->at(2)(kA, kB, v) * m_R->R(0, t,u,v);
                dVdQ(1) += m_Eab->at(0)(iA, iB, t)  * m_dEab_dQab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u,v);
                dVdQ(2) += m_Eab->at(0)(iA, iB, t)  * m_Eab->at(1)(jA, jB, u)  * m_dEab_dQab->at(2)(kA, kB, v) * m_R->R(0,t,u,v);
            }
        }
    }


    return 2. * M_PI / m_p * dVdQ * m_primitiveA->weight() * m_primitiveB->weight();
}

rowvec NuclearAttractionIntegralGD::PDerivative()
{
    rowvec dVdP = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();


    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){
                dVdP(0) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t+1,u,v);
                dVdP(1) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u+1,v);
                dVdP(2) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u,v+1);
            }
        }
    }

    return 2. * M_PI / m_p * dVdP * m_primitiveA->weight() * m_primitiveB->weight();
}



