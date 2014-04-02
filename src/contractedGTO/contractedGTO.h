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
    const vector<PrimitiveGTO>& primitiveGTOs() const;

    const rowvec &center() const;
    void setCenter(const rowvec* center);

private:
    vector<PrimitiveGTO> m_primitivesGTOs;
    const rowvec* m_center;
};
}

#endif // CONTRACTEDGTO_H
