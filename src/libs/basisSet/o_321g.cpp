#include "o_321g.h"

O_321G::O_321G()
{

    rowvec exp1 = {322.0370000, 48.430800, 10.4206000};
    rowvec exp2 = {7.4029400, 1.5762000};
    rowvec exp3 = {0.3736840};
    rowvec exp4 = {7.4029400, 1.5762000};
    rowvec exp5 = {0.3736840};

    rowvec weight1 = {0.0592394, 0.3515000, 0.7076580};
    rowvec weight2 = {-0.4044530, 1.2215600};
    rowvec weight3 = {1.0};
    rowvec weight4 = {0.2445860, 0.8539550};
    rowvec weight5 = {1.0};

    rowvec powers1 = {0, 0, 0};
    rowvec powers2 = {1, 0, 0};
    rowvec powers3 = {0, 1, 0};
    rowvec powers4 = {0, 0, 1};



    m_contractedGTO_1 = new ContractedGTO;
    m_contractedGTO_2 = new ContractedGTO;
    m_contractedGTO_3 = new ContractedGTO;
    m_contractedGTO_4 = new ContractedGTO;
    m_contractedGTO_5 = new ContractedGTO;
    m_contractedGTO_6 = new ContractedGTO;
    m_contractedGTO_7 = new ContractedGTO;
    m_contractedGTO_8 = new ContractedGTO;
    m_contractedGTO_9 = new ContractedGTO;


    //1s
    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exp1[0],weight1[0] ,powers1));
    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exp1[1],weight1[1] ,powers1));
    m_contractedGTO_1->addPrimitiveGTO(new PrimitiveGTO(exp1[2],weight1[2] ,powers1));


    //2s
    m_contractedGTO_2->addPrimitiveGTO(new PrimitiveGTO(exp2[0],weight2[0] ,powers1));
    m_contractedGTO_2->addPrimitiveGTO(new PrimitiveGTO(exp2[1],weight2[1] ,powers1));

    m_contractedGTO_3->addPrimitiveGTO(new PrimitiveGTO(exp3[0],weight3[0] ,powers1));

    //2px
    m_contractedGTO_4->addPrimitiveGTO(new PrimitiveGTO(exp4[0],weight4[0] ,powers2));
    m_contractedGTO_4->addPrimitiveGTO(new PrimitiveGTO(exp4[1],weight4[1] ,powers2));

    m_contractedGTO_5->addPrimitiveGTO(new PrimitiveGTO(exp5[0],weight5[0] ,powers2));

    //2py
    m_contractedGTO_6->addPrimitiveGTO(new PrimitiveGTO(exp4[0],weight4[0] ,powers3));
    m_contractedGTO_6->addPrimitiveGTO(new PrimitiveGTO(exp4[1],weight4[1] ,powers3));

    m_contractedGTO_7->addPrimitiveGTO(new PrimitiveGTO(exp5[0],weight5[0] ,powers3));

    //2pz
    m_contractedGTO_8->addPrimitiveGTO(new PrimitiveGTO(exp4[0],weight4[0] ,powers4));
    m_contractedGTO_8->addPrimitiveGTO(new PrimitiveGTO(exp4[1],weight4[1] ,powers4));

    m_contractedGTO_9->addPrimitiveGTO(new PrimitiveGTO(exp5[0],weight5[0] ,powers4));


    m_contractedGTOs.push_back(m_contractedGTO_1);
    m_contractedGTOs.push_back(m_contractedGTO_2);
    m_contractedGTOs.push_back(m_contractedGTO_3);
    m_contractedGTOs.push_back(m_contractedGTO_4);
    m_contractedGTOs.push_back(m_contractedGTO_5);
    m_contractedGTOs.push_back(m_contractedGTO_6);
    m_contractedGTOs.push_back(m_contractedGTO_7);
    m_contractedGTOs.push_back(m_contractedGTO_8);
    m_contractedGTOs.push_back(m_contractedGTO_9);


    m_angularMomentum = 1;

}

int O_321G::getAngularMomentum() const
{
    return m_angularMomentum;
}

