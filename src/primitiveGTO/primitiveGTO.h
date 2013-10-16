#ifndef PRIMITIVEGTO_H
#define PRIMITIVEGTO_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


class PrimitiveGTO
{
public:
    PrimitiveGTO(double alpha, double coff, string type);
    double alpha() const;
    double coff() const;

    void setAlpha(double alpha);
    void setCoff(double coff);

    string type() const;
    void setType(string type);

private:
    double m_alpha;
    double m_coff;
    string m_type;

};
#endif // PRIMITIVEGTO_H
