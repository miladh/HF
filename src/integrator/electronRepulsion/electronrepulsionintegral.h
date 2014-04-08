#ifndef ELECTRONREPULSIONINTEGRAL_H
#define ELECTRONREPULSIONINTEGRAL_H


#include <iostream>
#include <armadillo>
#include "../../math/hermitecoefficients.h"
#include "../../math/hermiteintegrals.h"
#include "../../primitiveGTO/primitiveGTO.h"


using namespace arma;
using namespace std;

namespace hf
{

class ElectronRepulsionIntegral
{
public:
    ElectronRepulsionIntegral(const int highestOrder,
                              const field<cube> *Eab, const field<cube> *Ecd,
                              const PrimitiveGTO *primitiveA,
                              const PrimitiveGTO *primitiveB, const PrimitiveGTO *primitiveC, const PrimitiveGTO *primitiveD);

    double evaluate();

private:
    HermiteIntegrals* m_R;
    const field<cube>* m_Eab;
    const field<cube>* m_Ecd;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
    const PrimitiveGTO* m_primitiveC;
    const PrimitiveGTO* m_primitiveD;

};
}
#endif // ELECTRONREPULSIONINTEGRAL_H
