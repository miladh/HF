#ifndef PRIMITIVEGTO_H
#define PRIMITIVEGTO_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


class PrimitiveGTO
{
public:
    PrimitiveGTO(const double &exponent = 0.0, const double &weight = 0.0 , const rowvec &powers = {0,0,0});

    const double &exponent() const;
    void setExponent(const double &exponent);

    const double &weight() const;
    void setWeight(const double &weight);

    const rowvec &powers() const;
    void setPowers(const rowvec &powers);

    const int &xPower() const;
    void setXPower(int xPower);

    const int &yPower() const;
    void setYPower(int yPower);

    const int &zPower() const;
    void setZPower(int zPower);

    const double &evaluate(const double &x, const double &y, const double &z);
private:
    double m_value;
    double m_exponent;
    double m_weight;
    int m_xPower, m_yPower, m_zPower;
    rowvec m_powers;


};
#endif // PRIMITIVEGTO_H
