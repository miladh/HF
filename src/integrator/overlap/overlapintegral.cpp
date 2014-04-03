#include "overlapintegral.h"

using namespace hf;

OverlapIntegral::OverlapIntegral(const field<cube>* Eab,
                                 const PrimitiveGTO* primitiveA,
                                 const PrimitiveGTO* primitiveB):
    m_Eab(Eab),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB)
{
//    field<cube>* A;
//    m_Eab = A;

//    cube AA;
//    m_Eab->at(0) = AA;

//    cout << "--------------"<< endl;
//    cout << "OVERLAP"  << endl;
//    cout << "m_Eab adr:   "  << m_Eab << endl;
//    cout << "m_cube1 adr:   " << &m_Eab->at(0) << endl;
//    cout << "m_cube2 adr:   " << &m_Eab->at(1) << endl;
//    cout << "m_cube3 adr:   " << &m_Eab->at(2) << endl;

//    cout << "--------------"<< endl;
//    cout << "OVERLAP"  << endl;
//    cout << "m_Eab adr:   "  << m_primitiveA << endl;
//    cout << "m_cube1 adr:   " << &m_primitiveA->center() << endl;

}


double OverlapIntegral::evaluate(int cor, int iA, int iB)
{
//    const cube& E =  m_Eab->at(cor);
    return m_Eab->at(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA->exponent() + m_primitiveB->exponent()));
}

double OverlapIntegral::evaluate()
{

    return    evaluate(0, m_primitiveA->xPower(), m_primitiveB->xPower())
            * evaluate(1, m_primitiveA->yPower(), m_primitiveB->yPower())
            * evaluate(2, m_primitiveA->zPower(), m_primitiveB->zPower())
            * m_primitiveA->weight() * m_primitiveB->weight();
}
