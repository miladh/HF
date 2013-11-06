#ifndef H_QuadZeta_H
#define H_QuadZeta_H

#include <basisSet/basisset.h>

class H_QuadZeta : public BasisSet
{
public:
    H_QuadZeta();

    int getAngularMomentum() const;


private:
    ContractedGTO m_contractedGTO_1;
    ContractedGTO m_contractedGTO_2;
    ContractedGTO m_contractedGTO_3;
    ContractedGTO m_contractedGTO_4;

};

#endif // H_QuadZeta_H
