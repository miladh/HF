#ifndef H_321G_H
#define H_321G_H
#include<basisSet/basisset.h>

class H_321G : public BasisSet
{
public:
    H_321G();
    int getAngularMomentum() const;

private:
    ContractedGTO m_contractedGTO_1;
    ContractedGTO m_contractedGTO_2;

};

#endif // H_321G_H
