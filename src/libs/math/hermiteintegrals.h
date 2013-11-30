#ifndef HERMITEINTEGRALS_H
#define HERMITEINTEGRALS_H

#include <iostream>
#include <armadillo>
#include<math/boys.h>


using namespace arma;
using namespace std;

class HermiteIntegrals
{
public:
    HermiteIntegrals(const int highestOrder);
    void setupR(const rowvec &PQ, const double &alpha, field<cube> &R, const rowvec &maxValues);

private:
    Boys *m_boys;
};

#endif // HERMITEINTEGRALS_H
