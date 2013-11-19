#ifndef C_321G_H
#define C_321G_H

#include<basisSet/splitValence/splitvalence.h>

class C_321G : public BasisSet
{
public:
    C_321G();

    int getAngularMomentum() const;

private:
    ContractedGTO m_contractedGTO_1;
    ContractedGTO m_contractedGTO_2;
    ContractedGTO m_contractedGTO_3;
    ContractedGTO m_contractedGTO_4;
    ContractedGTO m_contractedGTO_5;
    ContractedGTO m_contractedGTO_6;
    ContractedGTO m_contractedGTO_7;
    ContractedGTO m_contractedGTO_8;
    ContractedGTO m_contractedGTO_9;
};

#endif // C_321G_H
