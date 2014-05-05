#ifndef NUCLEARATTRACTIONINTEGRALGD_H
#define NUCLEARATTRACTIONINTEGRALGD_H



#include <iostream>
#include <armadillo>
#include "../../math/hermitecoefficients.h"
#include "../../math/hermiteintegrals.h"
#include "../../primitiveGTO/primitiveGTO.h"

using namespace arma;
using namespace std;

namespace hf
{

class NuclearAttractionIntegralGD
{
public:
    NuclearAttractionIntegralGD(const int highestOrder,
                                const field<cube> *Eab,
                                const field<cube> *dEab_dQab,
                                const PrimitiveGTO *primitiveA,
                                const PrimitiveGTO *primitiveB,
                                const rowvec* sourceCharge);

    void updateHermiteIntegrals();
    rowvec QDerivative();
    rowvec PDerivative();

private:
    HermiteIntegrals* m_R;
    const field<cube>* m_Eab;
    const field<cube>* m_dEab_dQab;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
    const rowvec* m_sourceCharge;
    double m_p;
    int m_tMax;
    int m_uMax;
    int m_vMax;

};

}

#endif // NUCLEARATTRACTIONINTEGRALGD_H
