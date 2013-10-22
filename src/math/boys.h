#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

class Boys
{
public:
    Boys();

    int factorial(int n);
    int doubleFactorial(int n);
    double boysFunction(double arg, int n);

private:
    double taylorExpandendBoys(double arg, int n, int nterms);
};

#endif // BOYSFUNCTION_H
