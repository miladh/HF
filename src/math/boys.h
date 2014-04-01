#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


namespace hf
{
class Boys
{
public:
    Boys(uint highestOrder);

    double getBoysFunctions(const int &n) const;
    rowvec getBoysFunctions() const;
    void evaluateBoysFunctions(const double &arg);

    static double factorial(const int &n);
    static double doubleFactorial(const int &n);

private:
    double m_arg;
    uint m_highestOrder;
    rowvec m_results;
    mat m_Ftabulated;

    void downwardRecursion();
    void readBoysForSmallArguments();
    double taylorExpandendBoys(uint nterms = 6, double dxt = 50.0/999) const;
    double asymptoticBoys() const;
};
}
#endif // BOYSFUNCTION_H
