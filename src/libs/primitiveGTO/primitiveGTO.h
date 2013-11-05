#ifndef PRIMITIVEGTO_H
#define PRIMITIVEGTO_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


class PrimitiveGTO
{
public:
    PrimitiveGTO(const double &exponent, const double &weight, const rowvec &powers);

    double exponent() const;
    void setExponent(const double &exponent);

    double weight() const;
    void setWeight(const double &weight);

    rowvec powers() const;
    void setPowers(const rowvec &powers);

private:
    double m_exponent;
    double m_weight;
    rowvec m_powers;

};
#endif // PRIMITIVEGTO_H
