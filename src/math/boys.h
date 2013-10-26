#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

class Boys
{
public:
    Boys(uint highestOrder);

    rowvec getBoysFunctions() const;
    void evaluateBoysFunctions(const double &arg);

    static double factorial(const int &n);
    static double doubleFactorial(const int &n);

private:
    double m_arg;
    uint m_highestOrder;
    rowvec m_results;
    mat m_Ftabulated;

    void results();
    void downwardRecursion();
    void readBoysForSmallArguments();
    double taylorExpandendBoys(uint nterms = 6, double dxt = 50.0/999) const;
    double asymptoticBoys() const;
};

#endif // BOYSFUNCTION_H
