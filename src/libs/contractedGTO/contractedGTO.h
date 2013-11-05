#ifndef CONTRACTEDGTO_H
#define CONTRACTEDGTO_H

#include <iostream>
#include <armadillo>
#include<primitiveGTO/primitiveGTO.h>

using namespace arma;
using namespace std;


class ContractedGTO
{
public:
    ContractedGTO();

    void addPrimitiveGTO(PrimitiveGTO *primitive);
    vector<PrimitiveGTO *> primitives() const;


    rowvec corePosition() const;
    void setCorePosition(const rowvec &corePosition);

    int getNumPrimitives() const;
private:
    vector<PrimitiveGTO *> m_primitives;
};


#endif // CONTRACTEDGTO_H
