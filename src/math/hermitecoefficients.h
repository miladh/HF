#ifndef HERMITECOEFFICIENTS_H
#define HERMITECOEFFICIENTS_H

#include <iostream>
#include <armadillo>

#include "../primitiveGTO/primitiveGTO.h"


using namespace arma;
using namespace std;

namespace hf
{
class HermiteCoefficients
{

public:
    HermiteCoefficients(const int maxAngularMomentum);
    HermiteCoefficients();

    const field<cube>* coefficients() const;
    const field<cube> *QDerivativeCoefficients() const;

    void updateE(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB, bool kin=true);
    void updatedE_dQ(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB, bool kin=true);

private:
    bool interiorPoint(int iA, int iB, int t);
    field<cube> m_E;
    field<cube> m_dE_dQ;

};
}
#endif // HERMITECOEFFICIENTS_H
