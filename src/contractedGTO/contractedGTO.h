#ifndef CONTRACTEDGTO_H
#define CONTRACTEDGTO_H

#include <iostream>
#include <armadillo>
#include "../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;

namespace hf
{
class ContractedGTO
{
public:
    ContractedGTO();

    void addPrimitive(PrimitiveGTO primitiveGTO);
    int getNumPrimitives() const;
    const PrimitiveGTO &getPrimitive(const int p) const;

    const rowvec &center() const;
    void setCenter(const rowvec &center);

    const vector<PrimitiveGTO>& primitivesGTOs() const;
    void setPrimitivesGTOs(const vector<PrimitiveGTO> &primitivesGTOs);

private:
    vector<PrimitiveGTO> m_primitivesGTOs;
    rowvec m_center;
};
}

#endif // CONTRACTEDGTO_H
