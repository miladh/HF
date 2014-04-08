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

    void updateR(const rowvec &PQ, const double &alpha);
    void updateR(const rowvec &PQ, const double &alpha,
                 const int tMax, const int uMax, const int vMax);

    double R(const int n, const int t, const int u, const int v) const;

private:
    Boys *m_boys;
    field<cube> m_R;
};
}
#endif // HERMITEINTEGRALS_H
