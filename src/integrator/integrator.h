#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>

#include<src/primitiveGTO/primitiveGTO.h>
#include<src/math/boys.h>

using namespace arma;
using namespace std;

class Integrator
{
public:
    Integrator();

    void setupE();
    void setupR(const rowvec &C);

    rowvec corePositionA() const;
    void setCorePositionA(const rowvec &corePositionA);

    rowvec corePositionB() const;
    void setCorePositionB(const rowvec &corePositionB);

    mat corePositions() const;
    void setCorePositions(const mat corePostions);

    uint maxAngularMomentum() const;
    void setMaxAngularMomentum(const uint &maxAngularMomentum);

    void addPrimitives(PrimitiveGTO *primitive);

    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);

    double boysFunction(double arg, int n);
    int factorial(int num);

private:
    rowvec m_corePositionA;
    rowvec m_corePositionB;
    mat m_corePositions;

    vector<PrimitiveGTO *> m_primitives;

    uint m_maxAngularMomentum;

    cube m_E[3]; // x,y,z cube
//    double ****m_R;

    vector<cube> m_R;




    bool interiorPoint(int iA, int iB, int t);

    double overlapIntegral(int cor, int iA, int iB);
    double kineticIntegral(int dim, int iA, int iB);

    double taylorExpandendBoys(double arg, int n, int nterms);
};

#endif // INTEGRATOR_H
