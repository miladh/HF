#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>

#include<math/hermitecoefficients.h>
#include<math/hermiteintegrals.h>

using namespace arma;
using namespace std;

class Integrator
{
public:
    Integrator();

    uint maxAngularMomentum() const;
    void setMaxAngularMomentum(const uint &maxAngularMomentum);

    double exponentA() const;
    void setExponentA(double exponentA);

    double exponentB() const;
    void setExponentB(double exponentB);

    double exponentC() const;
    void setExponentC(double exponentC);

    double exponentD() const;
    void setExponentD(double exponentD);


    rowvec corePositionA() const;
    void setCorePositionA(const rowvec &corePositionA);

    rowvec corePositionB() const;
    void setCorePositionB(const rowvec &corePositionB);

    rowvec corePositionC() const;
    void setCorePositionC(const rowvec &corePositionC);

    rowvec corePositionD() const;
    void setCorePositionD(const rowvec &corePositionD);


    double overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double nuclearAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB);
    double electronRepulsionIntegral(int iA, int jA, int kA, int iB, int jB, int kB,
                                     int iC, int jC, int kC, int iD, int jD, int kD);

    void updateHermiteCoefficients(bool oneParticleIntegral, bool twoParticleIntegral);
    void updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral);


   rowvec overlapIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
   rowvec kineticIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
   rowvec nuclearAttractionIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                               bool differentiateWrtA, bool differentiateWrtB, bool differentiateWrtC);


private:
    uint m_maxAngularMomentum;

    double m_exponentA;
    double m_exponentB;
    double m_exponentC;
    double m_exponentD;

    rowvec m_corePositionA;
    rowvec m_corePositionB;
    rowvec m_corePositionC;
    rowvec m_corePositionD;

    field<cube> m_Eab, m_Ecd;
    field<cube> m_dEab, m_dEcd;
    field<cube> m_Ree, m_Ren;

    HermiteCoefficients m_hermiteCoefficients;
    HermiteIntegrals *m_hermiteIntegrals;

    double overlapIntegral(int cor, int iA, int iB);
    double kineticIntegral(int cor, int iA, int iB);
    double overlapIntegral_derivative(int cor, int iA, int iB);
    double kineticIntegral_derivative(int cor, int iA, int iB);

    rowvec nuclearAttractionIntegral_R_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_P_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_C_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
};

#endif // INTEGRATOR_H

