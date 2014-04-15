#ifndef DIPOLEINTEGRAL_H
#define DIPOLEINTEGRAL_H

#include <iostream>
#include <armadillo>

#include "../../math/hermitecoefficients.h"
#include "../../primitiveGTO/primitiveGTO.h"

namespace hf
{

class DipoleIntegral
{
public:
    DipoleIntegral(const field<cube> *Eab,
                   const PrimitiveGTO *primitiveA,
                   const PrimitiveGTO *primitiveB);

    double evaluate(int cor, int iA, int iB, const double &P);
    double evaluate(int cor, int iA, int iB);
    rowvec evaluate();

    const PrimitiveGTO *primitiveA() const;
    const PrimitiveGTO *primitiveB() const;

private:
    const field<cube>* m_Eab;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
};
}
#endif // DIPOLEINTEGRAL_H
