#ifndef PRIMITIVEGTO_H
#define PRIMITIVEGTO_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


class PrimitiveGTO
{
public:
    PrimitiveGTO(double exponent, double weight);

    double exponent() const;
    void setExponent(double exponent);

    double weight() const;
    void setWeight(double weight);

private:
    double m_exponent;
    double m_weight;

};
#endif // PRIMITIVEGTO_H
