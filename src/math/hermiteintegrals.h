#ifndef HERMITEINTEGRALS_H
#define HERMITEINTEGRALS_H

#include <iostream>
#include <armadillo>
#include "../math/boys.h"


using namespace arma;
using namespace std;

namespace hf
{
class HermiteIntegrals
{
public:
    HermiteIntegrals(const int highestOrder);
    void setupR(const rowvec &PQ, const double &alpha, field<cube> &R);
    void setupR(const rowvec &PQ, const double &alpha, field<cube> &R,
                const int tMax, const int uMax, const int vMax);


private:
    Boys *m_boys;
};
}
#endif // HERMITEINTEGRALS_H
