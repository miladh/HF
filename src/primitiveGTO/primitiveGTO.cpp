#include "primitiveGTO.h"

PrimitiveGTO::PrimitiveGTO(double alpha, double coff, string type):
    m_alpha(alpha),
    m_coff(coff),
    m_type(type)
{
    setAlpha(m_alpha);
    setCoff(m_coff);
    setType(m_type);
}

double PrimitiveGTO::coff() const
{
    return m_coff;
}

void PrimitiveGTO::setCoff(double coff)
{
    m_coff = coff;
}


double PrimitiveGTO::alpha() const
{
    return m_alpha;
}

void PrimitiveGTO::setAlpha(double alpha)
{
    m_alpha = alpha;
}

string PrimitiveGTO::type() const
{
    return m_type;
}

void PrimitiveGTO::setType(string type)
{
    m_type = type;
}
