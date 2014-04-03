#include "overlapintegralgd.h"

using namespace hf;


OverlapIntegralGD::OverlapIntegralGD(OverlapIntegral *overlap,
                                     const field<cube> *dEab_dQab):
    m_overlap(overlap),
    m_dEab_dQab(dEab_dQab),
    m_primitiveA(overlap->primitiveA()),
    m_primitiveB(overlap->primitiveB())

{
//        cout << "--------------"<< endl;
//        cout << "OVERLAP"  << endl;
//        cout << "m_Eab adr:   "  << m_dEab_dQab << endl;
//        cout << "m_cube1 adr:   " << &m_dEab_dQab->at(0) << endl;
//        cout << "m_cube2 adr:   " << &m_dEab_dQab->at(1) << endl;
//        cout << "m_cube3 adr:   " << &m_dEab_dQab->at(2) << endl;

//        cout << "--------------"<< endl;
//        cout << "PrimOverlap"  << endl;
//        cout << "m_Eab adr:   "  << m_primitiveA << endl;
//        cout << "m_cube1 adr:   " << &m_primitiveA->center() << endl;
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


