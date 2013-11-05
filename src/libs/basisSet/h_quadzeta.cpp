#include "h_quadzeta.h"

H_QuadZeta::H_QuadZeta()
{
    rowvec exponent = {13.00773, 1.962079, 0.444529, 0.1219492};
    rowvec weight = {1.0, 1.0, 1.0, 1.0};
    rowvec powers = {0,0,0};


    m_contractedGTO_1 = new ContractedGTO;
    m_contractedGTO_2 = new ContractedGTO;
    m_contractedGTO_3 = new ContractedGTO;
    m_contractedGTO_4 = new ContractedGTO;

    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exponent[0],weight[0] ,powers));
    m_contractedGTO_2->addPrimitiveGTO(new PrimitiveGTO(exponent[1],weight[1] ,powers));
    m_contractedGTO_3->addPrimitiveGTO(new PrimitiveGTO(exponent[2],weight[2] ,powers));
    m_contractedGTO_4->addPrimitiveGTO(new PrimitiveGTO(exponent[3],weight[3] ,powers));

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
