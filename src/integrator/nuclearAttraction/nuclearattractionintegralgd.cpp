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

    double p = a + b;
    rowvec PC = (a*A + b*B)/p - C;

//    m_hermiteIntegrals->setupR(PC,p, m_Ren, iA+iB, jA+jB, kA+kB);
    m_R->updateR(PC,p);
}
rowvec NuclearAttractionIntegralGD::evaluate()
{

//    if(differentiateWrtA){
//        dVab += a/p * PDerivative(iA, jA, kA, iB, jB, kB)
//               + QDerivative(iA, jA, kA, iB, jB, kB) ;
//    }

//    if(differentiateWrtB){
//        dVab += b/p * PDerivative(iA, jA, kA, iB, jB, kB)
//                - QDerivative(iA, jA, kA, iB, jB, kB) ;
//    }

//    if(differentiateWrtC){
//        dVab -= CDerivative(iA, jA, kA, iB, jB, kB);
//    }

//    return 2 * M_PI / p * dVab;
}


rowvec NuclearAttractionIntegralGD::QDerivative()
{
    rowvec dVdQ = zeros<rowvec>(3);
    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    double p = a + b;

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdQ(0) += m_dEab_dQab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u)  * m_Eab->at(2)(kA, kB, v) * m_R->R(0, t,u,v);
                dVdQ(1) += m_Eab->at(0)(iA, iB, t)  * m_dEab_dQab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u,v);
                dVdQ(2) += m_Eab->at(0)(iA, iB, t)  * m_Eab->at(1)(jA, jB, u)  * m_dEab_dQab->at(2)(kA, kB, v) * m_R->R(0,t,u,v);
            }
        }
    }


    return 2 * M_PI / p * dVdQ;
}

rowvec NuclearAttractionIntegralGD::PDerivative()
{
    rowvec dVdP = zeros<rowvec>(3);
    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    double p = a + b;

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();


    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdP(0) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t+1,u,v);
                dVdP(1) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u+1,v);
                dVdP(2) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u,v+1);
            }
        }
    }

    return 2 * M_PI / p * dVdP;
}

rowvec NuclearAttractionIntegralGD::CDerivative()
{
    rowvec dVdC = zeros<rowvec>(3);
    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    double p = a + b;

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();


    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdC(0) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t+1,u,v);
                dVdC(1) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u+1,v);
                dVdC(2) += m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v) * m_R->R(0,t,u,v+1);
            }
        }
    }

    return 2 * M_PI / p * dVdC;
}


