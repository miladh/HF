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
    HermiteCoefficients();
    void setupE(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB, field<cube> &E, bool kin = true);

    void setup_dEdR(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB,
                    field<cube> &E, field<cube> &dE, bool kin=true);


private:
    bool interiorPoint(int iA, int iB, int t);

};
}
#endif // HERMITECOEFFICIENTS_H
