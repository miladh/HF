#include "h_321g.h"

H_321G::H_321G()
{

    rowvec exp1 = {5.4471780, 0.8245470};
    rowvec exp2 = {0.1831920};

    rowvec weight1 = {0.1562850, 0.9046910};
    rowvec weight2 = {1.0};

    rowvec powers = {0,0,0};


    m_contractedGTO_1 = new ContractedGTO;
    m_contractedGTO_2 = new ContractedGTO;


    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exp1[0],weight1[0] ,powers));
    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exp1[1],weight1[1] ,powers));
    m_contractedGTO_2->addPrimitiveGTO(new PrimitiveGTO(exp2[0],weight2[0] ,powers));

    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);

    m_angularMomentum = 0;

}

int H_321G::getAngularMomentum() const
{
    return m_angularMomentum;
}
