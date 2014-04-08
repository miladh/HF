#ifndef OVERLAPINTEGRALGD_H
#define OVERLAPINTEGRALGD_H


#include <iostream>
#include <armadillo>

#include "../../math/hermitecoefficients.h"
#include "../../primitiveGTO/primitiveGTO.h"
#include "overlapintegral.h"


namespace hf
{

class OverlapIntegralGD
{
public:
    OverlapIntegralGD(OverlapIntegral *overlap,
                      const field<cube>* dEab_dQab);

    double evaluate(int cor, int iA, int iB);
    rowvec evaluate();

    OverlapIntegral *overlap() const;

private:
    OverlapIntegral* m_overlap;
    const field<cube>* m_dEab_dQab;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
};

}
#endif // OVERLAPINTEGRALGD_H
