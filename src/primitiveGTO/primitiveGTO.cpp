#include "primitiveGTO.h"
using namespace hf;
PrimitiveGTO::PrimitiveGTO(const double &exponent, const double &weight, const rowvec &powers)
{
    setExponent(exponent);
    setWeight(weight);
    setPowers(powers);
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

const rowvec &PrimitiveGTO::center() const
{
    return (*m_center);
}

void PrimitiveGTO::setCenter(const rowvec *center)
{
    m_center = center;
}








