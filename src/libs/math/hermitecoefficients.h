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
    HermiteCoefficients(const double &a, const rowvec3 &A, const int &La,
                        const double &b, const rowvec3 &B, const int &Lb);


    inline field<cube> getCoefficients() const {return m_E;}

private:
    double m_a, m_b;
    rowvec3 m_A, m_B;
    int m_La, m_Lb;
    field<cube> m_E;

    void setupE();
    bool interiorPoint(int iA, int iB, int t);

};

#endif // HERMITECOEFFICIENTS_H
