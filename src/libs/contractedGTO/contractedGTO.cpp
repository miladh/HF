#include "contractedGTO.h"

ContractedGTO::ContractedGTO()
{
}

void ContractedGTO::addPrimitive(PrimitiveGTO primitiveGTO)
{
    m_primitivesGTOs.push_back(primitiveGTO);
}

const PrimitiveGTO &ContractedGTO::getPrimitive(const int p) const
{
    return m_primitivesGTOs.at(p);
}

int ContractedGTO::getNumPrimitives() const
{
    return m_primitivesGTOs.size();
}

double ContractedGTO::evaluate(const rowvec &R, const double &x, const double &y, const double &z) const
{
    double xDiff = x - R(0);
    double yDiff = y - R(1);
    double zDiff = z - R(2);
    double rSquared = (xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);
    double result = 0;
    for(uint i = 0; i < m_primitivesGTOs.size(); i++) {
        const PrimitiveGTO &p = m_primitivesGTOs.at(i);
        result += p.weight() * pow(xDiff, p.xPower()) *
                pow (yDiff, p.yPower()) * pow(zDiff, p.zPower())
                * exp(-p.exponent() * rSquared);
    }
    return result;
    return m_primitivesGTOs.size();
}
