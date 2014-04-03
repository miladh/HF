#ifndef NUCLEARATTRACTIONINTEGRAL_H
#define NUCLEARATTRACTIONINTEGRAL_H


#include <iostream>
#include <armadillo>
#include "../../math/hermitecoefficients.h"
#include "../../math/hermiteintegrals.h"
#include "../../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;

namespace hf
{

class NuclearAttractionIntegral
{
public:
    NuclearAttractionIntegral(HermiteIntegrals *hermiteIntegrals,
                              const field<cube> *Eab,
                              const PrimitiveGTO *primitiveA,
                              const PrimitiveGTO *primitiveB,
                              const rowvec* sourceCharge);

    double evaluate();
private:
    HermiteIntegrals* m_hermiteIntegrals;
    const field<cube>* m_Ren;
    const field<cube>* m_Eab;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
    const rowvec* m_sourceCharge;
};
}
#endif // NUCLEARATTRACTIONINTEGRAL_H
