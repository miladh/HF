#ifndef KINETICINTEGRALGD_H
#define KINETICINTEGRALGD_H

#include <iostream>
#include <armadillo>

#include "../../primitiveGTO/primitiveGTO.h"
#include "../overlap/overlapintegralgd.h"
#include "kineticintegral.h"

using namespace arma;
using namespace std;

namespace hf
{


class KineticIntegralGD
{
public:
    KineticIntegralGD(KineticIntegral* kinetic, OverlapIntegralGD *overlapGD);

    double evaluate(int cor, int iA, int iB);
    rowvec evaluate();

private:
    KineticIntegral* m_kinetic;
    OverlapIntegralGD* m_overlapGD;
    OverlapIntegral* m_overlap;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
};

}

#endif // KINETICINTEGRALGD_H
