#ifndef OVERLAPINTEGRAL_H
#define OVERLAPINTEGRAL_H

#include <iostream>
#include <armadillo>

#include "../../math/hermitecoefficients.h"
#include "../../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;

namespace hf
{

class OverlapIntegral
{
public:
    OverlapIntegral(const field<cube> *Eab,
                    const PrimitiveGTO *primitiveA,
                    const PrimitiveGTO *primitiveB);

    double evaluate(int cor, int iA, int iB);
    double evaluate();

    const PrimitiveGTO *primitiveA() const;
    const PrimitiveGTO *primitiveB() const;

private:
    const field<cube>* m_Eab;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
};

}
#endif // OVERLAPINTEGRAL_H
