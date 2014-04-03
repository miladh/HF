#include "nuclearattractionintegral.h"

using namespace hf;

NuclearAttractionIntegral::NuclearAttractionIntegral(const int highestOrder,
                                                     const field<cube> *Eab,
                                                     const PrimitiveGTO *primitiveA,
                                                     const PrimitiveGTO *primitiveB,
                                                     const rowvec* sourceCharge):
    m_R(new HermiteIntegrals(highestOrder)),
    m_Eab(Eab),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB),
    m_sourceCharge(sourceCharge)

{
//    cout << "--------------"<< endl;
//    cout << "NUC: "  << endl;
//    cout << "coreC adr   "  << m_sourceCharge<< endl;
//    cout << "primA adr   "  << m_primitiveB<< endl;
//    cout << "primB adr   "  << m_primitiveA<< endl;

//    cout << "--------------"<< endl;

//    cout <<"E adr:   " << &m_Eab << endl;
//    cout <<"E1 adr:   " << (m_Eab->at(0))(0,0,0) <<"    " << &m_Eab->at(0) << endl;
//    cout <<"E2 adr:   " << &(m_Eab->at(1)) << endl;
//    cout <<"E3 adr:   " << &(m_Eab->at(2)) << endl;

}

double NuclearAttractionIntegral::evaluate()
{
    const rowvec &A = m_primitiveA->center();
    const rowvec &B = m_primitiveB->center();
    const rowvec &C = (*m_sourceCharge);

    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();

    double p = a + b;
    rowvec PC = (a*A + b*B)/p - C;

    int tMax = m_primitiveA->xPower() + m_primitiveB->xPower() + 1;
    int uMax = m_primitiveA->yPower() + m_primitiveB->yPower() + 1;
    int vMax = m_primitiveA->zPower() + m_primitiveB->zPower() + 1;

    m_R->updateR(PC,p, tMax - 1 , uMax - 1, vMax - 1 );

    double result = 0.0;


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                result += m_Eab->at(0)(m_primitiveA->xPower(), m_primitiveB->xPower(), t)
                        * m_Eab->at(1)(m_primitiveA->yPower(), m_primitiveB->yPower(), u)
                        * m_Eab->at(2)(m_primitiveA->zPower(), m_primitiveB->zPower(), v)
                        * m_R->R(0,t,u,v);
            }
        }
    }

    return 2 * result * M_PI / p * m_primitiveA->weight() * m_primitiveB->weight();
}
