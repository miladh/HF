#include "primitiveGTO.h"

PrimitiveGTO::PrimitiveGTO(double exponent, double weight)
{
    setExponent(exponent);
    setWeight(weight);
}

double PrimitiveGTO::exponent() const
{
    return m_exponent;
}

void PrimitiveGTO::setExponent(double exponent)
{
    m_exponent = exponent;
}

double PrimitiveGTO::weight() const
{
    return m_weight;
}

void PrimitiveGTO::setWeight(double weight)
{
    m_weight = weight;
}


