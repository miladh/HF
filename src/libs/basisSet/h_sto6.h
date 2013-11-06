#ifndef H_STO6_H
#define H_STO6_H

#include<basisSet/basisset.h>
class H_STO6 : public BasisSet
{
public:
    H_STO6();
    int getAngularMomentum() const;

private:
    ContractedGTO *m_contractedGTO_1;
    ContractedGTO *m_contractedGTO_2;
    ContractedGTO *m_contractedGTO_3;
    ContractedGTO *m_contractedGTO_4;
    ContractedGTO *m_contractedGTO_5;
    ContractedGTO *m_contractedGTO_6;
};

#endif // H_STO6_H
