#ifndef KINETICINTEGRAL_H
#define KINETICINTEGRAL_H

#include <iostream>
#include <armadillo>

#include "../overlap/overlapintegral.h"
#include "../../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;

namespace hf
{

class KineticIntegral
{
public:
    KineticIntegral(OverlapIntegral *overlap,
                    const PrimitiveGTO *primitiveA,
                    const PrimitiveGTO *primitiveB);


    double evaluate(int cor, int iA, int iB);
    double evaluate();

private:
    OverlapIntegral* m_overlap;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
};

}
#endif // KINETICINTEGRAL_H
