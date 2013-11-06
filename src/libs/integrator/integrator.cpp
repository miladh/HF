#include "integrator.h"

Integrator::Integrator():
    m_exponentA(0),
    m_exponentB(0),
    m_exponentC(0),
    m_exponentD(0),
    m_corePositionA(rowvec(3)),
    m_corePositionB(rowvec(3)),
    m_corePositionC(rowvec(3)),
    m_corePositionD(rowvec(3))

{
        setMaxAngularMomentum(0);
}

void Integrator::setMaxAngularMomentum(const uint &maxAngularMomentum)
{
    m_maxAngularMomentum = maxAngularMomentum;
    m_Eab.set_size(3);
    m_Ecd.set_size(3);

    int nMax_en = 2 * m_maxAngularMomentum + 1;
    m_Ren.set_size(nMax_en);
    for(int n = 0; n < nMax_en; n++){
        m_Ren(n) = zeros(nMax_en, nMax_en, nMax_en);
    }

    int nMax_ee  = 4 * m_maxAngularMomentum + 1;
    m_Ree.set_size(nMax_ee);
    for(int n = 0; n < nMax_ee; n++){
        m_Ree(n) = zeros(nMax_ee, nMax_ee, nMax_ee);
    }

    m_hermiteIntegrals = new HermiteIntegrals(nMax_ee);
}

uint Integrator::maxAngularMomentum() const
{
    return m_maxAngularMomentum;
}


rowvec Integrator::corePositionA() const
{
    return m_corePositionA;
}

void Integrator::setCorePositionA(const rowvec &corePositionA)
{
    m_corePositionA = corePositionA;
}

rowvec Integrator::corePositionB() const
{
    return m_corePositionB;
}

void Integrator::setCorePositionB(const rowvec &corePositionB)
{
    m_corePositionB = corePositionB;
}


rowvec Integrator::corePositionC() const
{
    return m_corePositionC;
}

void Integrator::setCorePositionC(const rowvec &corePositionC)
{
    m_corePositionC = corePositionC;
}

rowvec Integrator::corePositionD() const
{
    return m_corePositionD;
}

void Integrator::setCorePositionD(const rowvec &corePositionD)
{
    m_corePositionD = corePositionD;
}


double Integrator::exponentA() const
{
    return m_exponentA;
}

void Integrator::setExponentA(double exponentA)
{
    m_exponentA = exponentA;
}
double Integrator::exponentB() const
{
    return m_exponentB;
}

void Integrator::setExponentB(double exponentB)
{
    m_exponentB = exponentB;
}
double Integrator::exponentC() const
{
    return m_exponentC;
}

void Integrator::setExponentC(double exponentC)
{
    m_exponentC = exponentC;
}
double Integrator::exponentD() const
{
    return m_exponentD;
}

void Integrator::setExponentD(double exponentD)
{
    m_exponentD = exponentD;
}


void Integrator::updateHermiteCoefficients(bool twoParticleIntegral)
{

    m_hermiteCoefficients.setupE(m_exponentA, m_corePositionA, m_maxAngularMomentum,
                                 m_exponentB, m_corePositionB, m_maxAngularMomentum,
                                 m_Eab);

    if(twoParticleIntegral){
        m_hermiteCoefficients.setupE(m_exponentC, m_corePositionC, m_maxAngularMomentum,
                                     m_exponentD, m_corePositionD, m_maxAngularMomentum,
                                     m_Ecd);
    }
}

/*---------------------------------------------------------------------------------------------------*/
double Integrator::overlapIntegral(int cor, int iA, int iB)
{
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    return m_Eab[cor](iA,iB,0) * sqrt(M_PI / p);
}


double Integrator::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {

    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}

/*---------------------------------------------------------------------------------------------------*/
double Integrator::kineticIntegral(int cor, int iA, int iB) {
    double b = m_exponentB;

    double S_iA_iBnn = overlapIntegral(cor, iA, iB + 2);
    double S_iA_iB = overlapIntegral(cor, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= overlapIntegral(cor, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB - 1) * S_iA_iBpp;
}

double Integrator::kineticIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {
    double T_iA_iB = kineticIntegral(0, iA, iB);
    double T_jA_jB = kineticIntegral(1, jA, jB);
    double T_kA_kB = kineticIntegral(2, kA, kB);

    double S_iA_iB = overlapIntegral(0, iA, iB);
    double S_jA_jB = overlapIntegral(1, jA, jB);
    double S_kA_kB = overlapIntegral(2, kA, kB);

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5;
    return result;
}

/*---------------------------------------------------------------------------------------------------*/
double Integrator::nuclearAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB)
{
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;

    double a = m_exponentA;
    double b = m_exponentB;

    double p = a + b;
    rowvec P  = (a*A + b*B)/p;
    rowvec PC = P - C;

    m_hermiteIntegrals->setupR(PC,p, m_Ren);

    double result = 0.0;

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
                result += m_Eab[0](iA, iB, t) * m_Eab[1](jA, jB, u) * m_Eab[2](kA, kB, v) * m_Ren(0)(t,u,v);
            }
        }
    }

    return 2 * result * M_PI / p;
}


/*---------------------------------------------------------------------------------------------------*/
double Integrator::electronRepulsionIntegral(int iA, int jA, int kA, int iB, int jB, int kB,
                                             int iC, int jC, int kC, int iD, int jD, int kD)
{
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;
    const rowvec &D = m_corePositionD;

    double a = m_exponentA;
    double b = m_exponentB;
    double c = m_exponentC;
    double d = m_exponentD;

    double p = a + b;
    double q = c + d;
    rowvec P  = (a*A + b*B)/p;
    rowvec Q  = (c*C + d*D)/q;

    double alpha = p*q/(p+q);
    rowvec PQ = P - Q;

    m_hermiteIntegrals->setupR(PQ,alpha, m_Ree);


    double result = 0.0;
    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;
    int kMax = iC + iD + 1;
    int lMax = jC + jD + 1;
    int mMax = kC + kD + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){

                double Etuv= m_Eab[0](iA, iB, t) * m_Eab[1](jA, jB, u) * m_Eab[2](kA, kB, v);

                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            double Eklm= m_Ecd[0](iC, iD, k) * m_Ecd[1](jC, jD, l) * m_Ecd[2](kC, kD, m);
                            result += Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                        }
                    }
                }


                result *=Etuv;
            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q));

    return result;

}







