#ifndef QUADZETA_H
#define QUADZETA_H

#include <basisSet/basisset.h>

class QuadZeta : public BasisSet
{
public:
    QuadZeta();

private:
    ContractedGTO *m_contractedGTO_1;
    ContractedGTO *m_contractedGTO_2;
    ContractedGTO *m_contractedGTO_3;
    ContractedGTO *m_contractedGTO_4;
};

#endif // QUADZETA_H
