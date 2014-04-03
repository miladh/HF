#ifndef ELECTRONREPULSIONINTEGRALGD_H
#define ELECTRONREPULSIONINTEGRALGD_H


#include <iostream>
#include <armadillo>
#include "../../math/hermitecoefficients.h"
#include "../../math/hermiteintegrals.h"
#include "../../primitiveGTO/primitiveGTO.h"


using namespace arma;
using namespace std;

namespace hf
{


class ElectronRepulsionIntegralGD
{
public:
    ElectronRepulsionIntegralGD(const int highestOrder,
                                const field<cube> *Eab,
                                const field<cube> *Ecd,
                                const PrimitiveGTO *primitiveA,
                                const PrimitiveGTO *primitiveB,
                                const PrimitiveGTO *primitiveC,
                                const PrimitiveGTO *primitiveD);


    void updateHermiteIntegrals();
    rowvec QabDerivative();
    rowvec PabDerivative();
    rowvec QcdDerivative();
    rowvec PcdDerivative();
private:
    HermiteIntegrals* m_R;
    const field<cube>* m_Eab;
    const field<cube>* m_Ecd;
    const field<cube>* m_dEab_dQab;
    const field<cube>* m_dEcd_dQcd;
    const PrimitiveGTO* m_primitiveA;
    const PrimitiveGTO* m_primitiveB;
    const PrimitiveGTO* m_primitiveC;
    const PrimitiveGTO* m_primitiveD;

};

}
#endif // ELECTRONREPULSIONINTEGRALGD_H
