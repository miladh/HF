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
    void setupE(const double &a, const double &b,
                const rowvec3 &R, field<cube> &E, bool kin = true);

    void setup_dEdR(const double &a, const double &b, const rowvec3 &R,
                    field<cube> &E, field<cube> &dE, bool kin=true);


private:
    bool interiorPoint(int iA, int iB, int t);

};

#endif // HERMITECOEFFICIENTS_H
