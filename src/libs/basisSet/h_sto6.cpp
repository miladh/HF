#include "h_sto6.h"

H_STO6::H_STO6()
{

    rowvec exponent = {35.52322122, 6.513143725, 1.822142904,
                       0.625955266, 0.243076747, 0.100112428};
    rowvec weight = {0.00916359628,0.04936149294,0.16853830490,
                    0.37056279970, 0.41649152980,0.13033408410} ;

    rowvec powers = {0,0,0};


    m_contractedGTO_1 = new ContractedGTO;
    m_contractedGTO_2 = new ContractedGTO;
    m_contractedGTO_3 = new ContractedGTO;
    m_contractedGTO_4 = new ContractedGTO;
    m_contractedGTO_5 = new ContractedGTO;
    m_contractedGTO_6 = new ContractedGTO;


    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exponent[0],weight[0] ,powers));
    m_contractedGTO_2->addPrimitiveGTO(new PrimitiveGTO(exponent[1],weight[1] ,powers));
    m_contractedGTO_3->addPrimitiveGTO(new PrimitiveGTO(exponent[2],weight[2] ,powers));
    m_contractedGTO_4->addPrimitiveGTO(new PrimitiveGTO(exponent[3],weight[3] ,powers));
    m_contractedGTO_5->addPrimitiveGTO(new PrimitiveGTO(exponent[4],weight[4] ,powers));
    m_contractedGTO_6->addPrimitiveGTO(new PrimitiveGTO(exponent[5],weight[5] ,powers));



    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);
    m_contractedGTOs.push_back(m_contractedGTO_3);
    m_contractedGTOs.push_back(m_contractedGTO_4);
    m_contractedGTOs.push_back(m_contractedGTO_5);
    m_contractedGTOs.push_back(m_contractedGTO_6);

    m_angularMomentum = 0;
}


int H_STO6::getAngularMomentum() const
{
    return m_angularMomentum;
}
