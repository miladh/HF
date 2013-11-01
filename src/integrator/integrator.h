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
    void setupR(const rowvec &PQ, const double &alpha, field<cube> &R);

    rowvec corePositionA() const;
    void setCorePositionA(const rowvec &corePositionA);

    rowvec corePositionB() const;
    void setCorePositionB(const rowvec &corePositionB);

    rowvec corePositionC() const;
    void setCorePositionC(const rowvec &corePositionC);

    rowvec corePositionD() const;
    void setCorePositionD(const rowvec &corePositionD);

    uint maxAngularMomentum() const;
    void setMaxAngularMomentum(const uint &maxAngularMomentum);

    void addPrimitives(PrimitiveGTO *primitive);

    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);  
    double nuclearAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double electronRepulsionIntegral(int iA, int jA, int kA, int iB, int jB, int kB,
                                     int iC, int jC, int kC, int iD, int jD, int kD);

    double exponentA() const;
    void setExponentA(double exponentA);

    double exponentB() const;
    void setExponentB(double exponentB);

    double exponentC() const;
    void setExponentC(double exponentC);

    double exponentD() const;
    void setExponentD(double exponentD);

private:
    rowvec m_corePositionA;
    rowvec m_corePositionB;
    rowvec m_corePositionC;
    rowvec m_corePositionD;

    double m_exponentA;
    double m_exponentB;
    double m_exponentC;
    double m_exponentD;

    vector<PrimitiveGTO *> m_primitives;

    uint m_maxAngularMomentum;

    cube m_E[3]; // x,y,z cube
    vector<cube> m_R;
    field<cube> m_Ree;
    field<cube> m_Ren;


    bool interiorPoint(int iA, int iB, int t);

    double overlapIntegral(int cor, int iA, int iB);
    double kineticIntegral(int cor, int iA, int iB);

    Boys* m_boys;
};

#endif // INTEGRATOR_H
