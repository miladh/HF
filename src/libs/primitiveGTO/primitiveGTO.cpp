#include "primitiveGTO.h"

PrimitiveGTO::PrimitiveGTO(const double &exponent, const double &weight, const rowvec &powers)
{
    setExponent(exponent);
    setWeight(weight);
    setPowers(powers);
}


const double& PrimitiveGTO::evaluate(const double &x,const double &y,const double &z)
{
    m_value = m_weight * pow(x, m_xPower) * pow (y, m_yPower) * pow(z, m_zPower)
            * exp(-m_exponent * (x*x + y*y + z*z));

    return m_value;
}

const double& PrimitiveGTO::exponent() const
{
    return m_exponent;
}

void PrimitiveGTO::setExponent(const double &exponent)
{
    m_exponent = exponent;
}

const double& PrimitiveGTO::weight() const
{
    return m_weight;
}

void PrimitiveGTO::setWeight(const double &weight)
{
    m_weight = weight;
}

const rowvec& PrimitiveGTO::powers() const
{
    return m_powers;
}

void PrimitiveGTO::setPowers(const rowvec &powers)
{
    m_powers = powers;
    setXPower(m_powers(0));
    setYPower(m_powers(1));
    setZPower(m_powers(2));
}


const int& PrimitiveGTO::xPower() const
{
    return m_xPower;
}
void PrimitiveGTO::setXPower(int xPower)
{
    m_xPower = xPower;
}


const int& PrimitiveGTO::yPower() const
{
    return m_yPower;
}
void PrimitiveGTO::setYPower(int yPower)
{
    m_yPower = yPower;
}


const int& PrimitiveGTO::zPower() const
{
    return m_zPower;
}

void PrimitiveGTO::setZPower(int zPower)
{
    m_zPower = zPower;
}







