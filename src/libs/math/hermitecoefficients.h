#ifndef HERMITECOEFFICIENTS_H
#define HERMITECOEFFICIENTS_H

#include <iostream>
#include <armadillo>


using namespace arma;
using namespace std;

class HermiteCoefficients
{

public:
    HermiteCoefficients();
    void setupE(const double &a, const rowvec3 &A, const int &La,
                const double &b, const rowvec3 &B, const int &Lb, field<cube> &E);

private:
    bool interiorPoint(int iA, int iB, int t);

};

#endif // HERMITECOEFFICIENTS_H
