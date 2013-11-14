#include "h_quadzeta.h"

H_QuadZeta::H_QuadZeta()
{
//    rowvec exponent = {13.00773, 1.962079, 0.444529, 0.1219492};
    rowvec exponent = {2.962079, 1.962079, 1.962079, 1.962079};
    rowvec weight = {1.0, 1.0, 1.0, 1.0};
    rowvec powers = {0,0,0};



    PrimitiveGTO primitiveGTO_1(exponent(0),weight(0) ,powers);
    PrimitiveGTO primitiveGTO_2(exponent(1),weight(1) ,powers);
    PrimitiveGTO primitiveGTO_3(exponent(2),weight(2) ,powers);
    PrimitiveGTO primitiveGTO_4(exponent(3),weight(3) ,powers);

    m_contractedGTO_1.addPrimitive(primitiveGTO_1);
    m_contractedGTO_2.addPrimitive(primitiveGTO_2);
    m_contractedGTO_3.addPrimitive(primitiveGTO_3);
    m_contractedGTO_4.addPrimitive(primitiveGTO_4);

    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);
    m_contractedGTOs.push_back(m_contractedGTO_3);
    m_contractedGTOs.push_back(m_contractedGTO_4);

    m_angularMomentum = 0;

}


int H_QuadZeta::getAngularMomentum() const
{
    return m_angularMomentum;
}
