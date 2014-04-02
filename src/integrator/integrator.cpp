#include "integrator.h"

using namespace hf;

Integrator::Integrator():
    m_corePositionC(rowvec(3))
{
}

void Integrator::setMaxAngularMomentum(const int maxAngularMomentum)
{
    m_Eab.set_size(3);  m_Ecd.set_size(3);
    m_dEab.set_size(3); m_dEcd.set_size(3);
    int iAmax = maxAngularMomentum + 3;
    int iBmax = maxAngularMomentum + 3;
    int tmax  = iAmax + iBmax - 1;
    for(int cor = 0; cor < 3; cor++){
        m_Eab(cor)  = zeros(iAmax, iBmax, tmax);
        m_Ecd(cor)  = zeros(iAmax, iBmax, tmax);
        m_dEab(cor) = zeros(iAmax, iBmax, tmax);
        m_dEcd(cor) = zeros(iAmax, iBmax, tmax);
    }

    int nMax_en = 2 * maxAngularMomentum + 1;
    m_Ren.set_size(nMax_en);
    for(int n = 0; n < nMax_en; n++){
        m_Ren(n) = zeros(nMax_en, nMax_en, nMax_en);
    }

    int nMax_ee  = 4 * maxAngularMomentum + 1;
    m_Ree.set_size(nMax_ee);
    for(int n = 0; n < nMax_ee; n++){
        m_Ree(n) = zeros(nMax_ee, nMax_ee, nMax_ee);
    }

    m_hermiteIntegrals = new HermiteIntegrals(nMax_ee);


}

//rowvec Integrator::corePositionC() const
//{
//    return m_primitiveC.center();
//}

void Integrator::setCorePositionC(const rowvec &corePositionC)
{
    m_corePositionC = corePositionC;
}

void Integrator::setPrimitiveA(const PrimitiveGTO &primitiveA)
{
    m_primitiveA = primitiveA;
}

void Integrator::setPrimitiveB(const PrimitiveGTO &primitiveB)
{
    m_primitiveB = primitiveB;
}

void Integrator::setPrimitiveC(const PrimitiveGTO &primitiveC)
{
    m_primitiveC = primitiveC;
}

void Integrator::setPrimitiveD(const PrimitiveGTO &primitiveD)
{
    m_primitiveD = primitiveD;
}

void Integrator::updateHermiteCoefficients(bool oneParticleIntegral, bool twoParticleIntegral, bool kin)
{

    if(oneParticleIntegral){
        if(kin){
            m_hermiteCoefficients.setupE(m_primitiveA, m_primitiveB, m_Eab);
        }else{
            m_hermiteCoefficients.setupE(m_primitiveA, m_primitiveB, m_Eab,false);
        }

    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setupE(m_primitiveC, m_primitiveD, m_Ecd, false);
    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }
}

void Integrator::updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral,bool kin)
{

    if(oneParticleIntegral){
        if(kin){
        m_hermiteCoefficients.setup_dEdR(m_primitiveA, m_primitiveB, m_Eab, m_dEab);
        }else{
            m_hermiteCoefficients.setup_dEdR(m_primitiveA, m_primitiveB, m_Eab, m_dEab, false);
        }

    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setup_dEdR(m_primitiveC, m_primitiveD, m_Ecd,m_dEcd, false);
    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }

}

/********************************************************************************************
 *
 *                                  Molecular Gaussian Integrals
 *
 * ******************************************************************************************/

double Integrator::overlapIntegral(int cor, int iA, int iB)
{
    return m_Eab(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA.exponent() + m_primitiveB.exponent()));
}

double Integrator::overlapIntegral()
{
    return    overlapIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower())
            * overlapIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower())
            * overlapIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower())
            * m_primitiveA.weight() * m_primitiveB.weight();
}

/*---------------------------------------------------------------------------------------------------*/

double Integrator::kineticIntegral(int cor, int iA, int iB) {
    double b = m_primitiveB.exponent();

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

double Integrator::kineticIntegral() {
    double T_iA_iB = kineticIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double T_jA_jB = kineticIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double T_kA_kB = kineticIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    double S_iA_iB = overlapIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double S_jA_jB = overlapIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double S_kA_kB = overlapIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    double result = T_iA_iB * S_jA_jB * S_kA_kB + S_iA_iB * T_jA_jB * S_kA_kB + S_iA_iB * S_jA_jB * T_kA_kB;
    result *= -0.5 * m_primitiveA.weight() * m_primitiveB.weight();
    return result;
}

/*---------------------------------------------------------------------------------------------------*/

double Integrator::nuclearAttractionIntegral()
{
    const rowvec &A = m_primitiveA.center();
    const rowvec &B = m_primitiveB.center();
    const rowvec &C = m_corePositionC;

    const double &a  = m_primitiveA.exponent();
    const double &b  = m_primitiveB.exponent();

    double p = a + b;
    rowvec PC = (a*A + b*B)/p - C;

    int tMax = m_primitiveA.xPower() + m_primitiveB.xPower() + 1;
    int uMax = m_primitiveA.yPower() + m_primitiveB.yPower() + 1;
    int vMax = m_primitiveA.zPower() + m_primitiveB.zPower() + 1;

    m_hermiteIntegrals->setupR(PC,p, m_Ren, tMax - 1 , uMax - 1, vMax - 1 );

    double result = 0.0;


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                result += m_Eab(0)(m_primitiveA.xPower(), m_primitiveB.xPower(), t)
                        * m_Eab(1)(m_primitiveA.yPower(), m_primitiveB.yPower(), u)
                        * m_Eab(2)(m_primitiveA.zPower(), m_primitiveB.zPower(), v)
                        * m_Ren(0)(t,u,v);
            }
        }
    }

    return 2 * result * M_PI / p * m_primitiveA.weight() * m_primitiveB.weight();
}

/*---------------------------------------------------------------------------------------------------*/

double Integrator::electronRepulsionIntegral()
{
    const rowvec &A = m_primitiveA.center();
    const rowvec &B = m_primitiveB.center();
    const rowvec &C = m_primitiveC.center();
    const rowvec &D = m_primitiveD.center();

    const double &a  = m_primitiveA.exponent();
    const double &b  = m_primitiveB.exponent();
    const double &c  = m_primitiveC.exponent();
    const double &d  = m_primitiveD.exponent();

    double p = a + b;
    double q = c + d;

    double alpha = p*q/(p+q);
    rowvec PQ = (a*A + b*B)/p - (c*C + d*D)/q;


    double result = 0.0;
    int tMax = m_primitiveA.xPower() + m_primitiveB.xPower() + 1;
    int uMax = m_primitiveA.yPower() + m_primitiveB.yPower() + 1;
    int vMax = m_primitiveA.zPower() + m_primitiveB.zPower() + 1;
    int kMax = m_primitiveC.xPower() + m_primitiveD.xPower() + 1;
    int lMax = m_primitiveC.yPower() + m_primitiveD.yPower() + 1;
    int mMax = m_primitiveC.zPower() + m_primitiveD.zPower() + 1;

    m_hermiteIntegrals->setupR(PQ,alpha, m_Ree, tMax + kMax - 2,
                               uMax + lMax - 2, vMax + mMax - 2);


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){

                double Etuv = m_Eab(0)(m_primitiveA.xPower() , m_primitiveB.xPower(), t)
                            * m_Eab(1)(m_primitiveA.yPower() , m_primitiveB.yPower(), u)
                            * m_Eab(2)(m_primitiveA.zPower() , m_primitiveB.zPower(), v);

                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            double Eklm = m_Ecd(0)(m_primitiveC.xPower() , m_primitiveD.xPower(), k)
                                        * m_Ecd(1)(m_primitiveC.yPower() , m_primitiveD.yPower(), l)
                                        * m_Ecd(2)(m_primitiveC.zPower() , m_primitiveD.zPower(), m);
                            result += Etuv * Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                        }
                    }
                }

            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q))
            * m_primitiveA.weight() * m_primitiveB.weight()
            * m_primitiveC.weight() * m_primitiveD.weight();

    return result;

}


/********************************************************************************************
 *
 *                  Molecular Gaussian Integral Geometrical Derivatives (GD)
 *
 * ******************************************************************************************/

double Integrator::overlapIntegral_derivative(int cor, int iA, int iB)
{
    return m_dEab(cor)(iA,iB,0) * sqrt(M_PI / (m_primitiveA.exponent() + m_primitiveB.exponent()));
}


rowvec Integrator::overlapIntegral_derivative()
{
    rowvec dS = zeros<rowvec>(3);

    dS(0) =   overlapIntegral_derivative(0, m_primitiveA.xPower(), m_primitiveB.xPower())
            * overlapIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower())
            * overlapIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    dS(1) =   overlapIntegral_derivative(1, m_primitiveA.yPower(), m_primitiveB.yPower())
            * overlapIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower())
            * overlapIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    dS(2) =   overlapIntegral_derivative(2, m_primitiveA.zPower(), m_primitiveB.zPower())
            * overlapIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower())
            * overlapIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower());

    return dS;
}

/*---------------------------------------------------------------------------------------------------*/

double Integrator::kineticIntegral_derivative(int cor, int iA, int iB) {
    double b = m_primitiveB.exponent();

    double dS_iA_iBnn = overlapIntegral_derivative(cor, iA, iB + 2);
    double dS_iA_iB   = overlapIntegral_derivative(cor, iA, iB);
    double dS_iA_iBpp;

    if(iB - 2 >= 0) {
        dS_iA_iBpp= overlapIntegral_derivative(cor, iA, iB - 2);
    } else {
        dS_iA_iBpp = 0;
    }
    return 4 * b * b * dS_iA_iBnn - 2*b * (2*iB + 1) * dS_iA_iB + iB * (iB - 1) * dS_iA_iBpp;
}


rowvec Integrator::kineticIntegral_derivative() {
    rowvec dT = zeros<rowvec>(3);

    double T_iA_iB = kineticIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double T_jA_jB = kineticIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double T_kA_kB = kineticIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    double dT_iA_iB = kineticIntegral_derivative(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double dT_jA_jB = kineticIntegral_derivative(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double dT_kA_kB = kineticIntegral_derivative(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    double S_iA_iB = overlapIntegral(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double S_jA_jB = overlapIntegral(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double S_kA_kB = overlapIntegral(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    double dS_iA_iB = overlapIntegral_derivative(0, m_primitiveA.xPower(), m_primitiveB.xPower());
    double dS_jA_jB = overlapIntegral_derivative(1, m_primitiveA.yPower(), m_primitiveB.yPower());
    double dS_kA_kB = overlapIntegral_derivative(2, m_primitiveA.zPower(), m_primitiveB.zPower());

    dT(0) = (dT_iA_iB * S_jA_jB * S_kA_kB) + (dS_iA_iB * T_jA_jB * S_kA_kB) + (dS_iA_iB * S_jA_jB  * T_kA_kB);
    dT(1) = (T_iA_iB * dS_jA_jB * S_kA_kB) + (S_iA_iB * dT_jA_jB * S_kA_kB) + (S_iA_iB * dS_jA_jB  * T_kA_kB);
    dT(2) = (T_iA_iB * S_jA_jB * dS_kA_kB) + (S_iA_iB * T_jA_jB * dS_kA_kB) + (S_iA_iB * S_jA_jB  * dT_kA_kB);

    dT *= -0.5;
    return dT;
}

/*---------------------------------------------------------------------------------------------------*/

rowvec Integrator::nuclearAttractionIntegral_R_derivative(int iA, int jA, int kA, int iB, int jB, int kB)
{
    rowvec dVdR = zeros<rowvec>(3);

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdR(0) += m_dEab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u)  * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v);
                dVdR(1) += m_Eab(0)(iA, iB, t)  * m_dEab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v);
                dVdR(2) += m_Eab(0)(iA, iB, t)  * m_Eab(1)(jA, jB, u)  * m_dEab(2)(kA, kB, v) * m_Ren(0)(t,u,v);
            }
        }
    }


    return dVdR;
}


rowvec Integrator::nuclearAttractionIntegral_P_derivative(int iA, int jA, int kA, int iB, int jB, int kB)
{
    rowvec dVdP = zeros<rowvec>(3);

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdP(0) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t+1,u,v);
                dVdP(1) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u+1,v);
                dVdP(2) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v+1);
            }
        }
    }

    return dVdP;
}


rowvec Integrator::nuclearAttractionIntegral_C_derivative(int iA, int jA, int kA, int iB, int jB, int kB)
{
    rowvec dVdC = zeros<rowvec>(3);

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;


    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){
                dVdC(0) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t+1,u,v);
                dVdC(1) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u+1,v);
                dVdC(2) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v+1);
            }
        }
    }

    return dVdC;
}

rowvec Integrator::nuclearAttractionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                                        bool differentiateWrtC)
{
    rowvec dVab = zeros<rowvec>(3);
    const rowvec &A = m_primitiveA.center();
    const rowvec &B = m_primitiveB.center();
    const rowvec &C = m_corePositionC;

    const double &a  = m_primitiveA.exponent();
    const double &b  = m_primitiveB.exponent();

    double p = a + b;
    rowvec PC = (a*A + b*B)/p - C;


    int iA = m_primitiveA.xPower();
    int jA = m_primitiveA.yPower();
    int kA = m_primitiveA.zPower();
    int iB = m_primitiveB.xPower();
    int jB = m_primitiveB.yPower();
    int kB = m_primitiveB.zPower();

//    m_hermiteIntegrals->setupR(PC,p, m_Ren, iA+iB, jA+jB, kA+kB);
    m_hermiteIntegrals->setupR(PC,p, m_Ren);

    if(differentiateWrtA){
        dVab += a/p * nuclearAttractionIntegral_P_derivative(iA, jA, kA, iB, jB, kB)
               + nuclearAttractionIntegral_R_derivative(iA, jA, kA, iB, jB, kB) ;
    }

    if(differentiateWrtB){
        dVab += b/p * nuclearAttractionIntegral_P_derivative(iA, jA, kA, iB, jB, kB)
                - nuclearAttractionIntegral_R_derivative(iA, jA, kA, iB, jB, kB) ;
    }

    if(differentiateWrtC){
        dVab -= nuclearAttractionIntegral_C_derivative(iA, jA, kA, iB, jB, kB);
    }

        return 2 * M_PI / p * dVab;
}

/*---------------------------------------------------------------------------------------------------*/

rowvec Integrator::electronRepulsionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                                        bool differentiateWrtC, bool differentiateWrtD)
{
    rowvec dQabcd = zeros<rowvec>(3);
    rowvec dQdPab, dQdPcd;
    rowvec dQdRab, dQdRcd;
    const rowvec &A = m_primitiveA.center();
    const rowvec &B = m_primitiveB.center();
    const rowvec &C = m_primitiveC.center();
    const rowvec &D = m_primitiveD.center();

    const double &a  = m_primitiveA.exponent();
    const double &b  = m_primitiveB.exponent();
    const double &c  = m_primitiveC.exponent();
    const double &d  = m_primitiveD.exponent();

    double p = a + b;
    double q = c + d;
    double alpha = p*q/(p+q);
    rowvec PQ = (a*A + b*B)/p - (c*C + d*D)/q;


    int iA = m_primitiveA.xPower();
    int jA = m_primitiveA.yPower();
    int kA = m_primitiveA.zPower();
    int iB = m_primitiveB.xPower();
    int jB = m_primitiveB.yPower();
    int kB = m_primitiveB.zPower();
    int iC = m_primitiveC.xPower();
    int jC = m_primitiveC.yPower();
    int kC = m_primitiveC.zPower();
    int iD = m_primitiveD.xPower();
    int jD = m_primitiveD.yPower();
    int kD = m_primitiveD.zPower();


//    m_hermiteIntegrals->setupR(PQ,alpha, m_Ree, iA+iB+iC+iD,
//                               jA+jB+jC+jD, kA+kB+kC+kD);
    m_hermiteIntegrals->setupR(PQ,alpha, m_Ree);


    if(differentiateWrtA || differentiateWrtB){
        dQdPab = electronRepulsionIntegral_Pab_derivative(iA, jA, kA, iB, jB, kB,
                                                      iC, jC, kC, iD, jD, kD);
        dQdRab = electronRepulsionIntegral_Rab_derivative(iA, jA, kA, iB, jB, kB,
                                                          iC, jC, kC, iD, jD, kD) ;
    }
    if(differentiateWrtC || differentiateWrtD){
        dQdRcd = electronRepulsionIntegral_Rcd_derivative(iA, jA, kA, iB, jB, kB,
                                                          iC, jC, kC, iD, jD, kD);
        dQdPcd = electronRepulsionIntegral_Pcd_derivative(iA, jA, kA, iB, jB, kB,
                                                          iC, jC, kC, iD, jD, kD);
    }


    if(differentiateWrtA){
        dQabcd  += a/p * dQdPab + dQdRab;
    }

    if(differentiateWrtB){
        dQabcd  += b/p * dQdPab - dQdRab;
    }

    if(differentiateWrtC){
        dQabcd  += c/q * dQdPcd + dQdRcd;
    }

    if(differentiateWrtD){
        dQabcd  += d/q * dQdPcd - dQdRcd;
    }

    return dQabcd * 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q));

}

rowvec Integrator::electronRepulsionIntegral_Pab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        int iC, int jC, int kC, int iD, int jD, int kD)
{
    rowvec dQdPab = zeros<rowvec>(3);

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

                            dQdPab(0) += Etuv*Eklm  * m_Ree(0)(t+1+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dQdPab(1) += Etuv*Eklm  * m_Ree(0)(t+k,u+1+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dQdPab(2) += Etuv*Eklm  * m_Ree(0)(t+k,u+l,v+1+m) * (1 - 2* ((k+l+m)%2));

                        }
                    }
                }

            }
        }
    }

    return dQdPab;

}


rowvec Integrator::electronRepulsionIntegral_Pcd_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        int iC, int jC, int kC, int iD, int jD, int kD)
{
    rowvec dQdPcd = zeros<rowvec>(3);

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

                            dQdPcd(0) += Etuv*Eklm  * m_Ree(0)(t+k+1,u+l,v+m) * (1 - 2* ((k+1+l+m)%2));
                            dQdPcd(1) += Etuv*Eklm  * m_Ree(0)(t+k,u+l+1,v+m) * (1 - 2* ((k+l+1+m)%2));
                            dQdPcd(2) += Etuv*Eklm  * m_Ree(0)(t+k,u+l,v+m+1) * (1 - 2* ((k+l+m+1)%2));

                        }
                    }
                }

            }
        }
    }

    return dQdPcd;

}

rowvec Integrator::electronRepulsionIntegral_Rab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        int iC, int jC, int kC, int iD, int jD, int kD)
{

    rowvec dQdRab = zeros<rowvec>(3);

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;
    int kMax = iC + iD + 1;
    int lMax = jC + jD + 1;
    int mMax = kC + kD + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){

                rowvec3 dEtuv = zeros<rowvec>(3);;
                dEtuv(0) = m_dEab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v);
                dEtuv(1) = m_Eab(0)(iA, iB, t) * m_dEab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v);
                dEtuv(2) = m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_dEab(2)(kA, kB, v);

                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            double Eklm= m_Ecd(0)(iC, iD, k) * m_Ecd(1)(jC, jD, l) * m_Ecd(2)(kC, kD, m);

                            dQdRab(0) += dEtuv(0) * Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dQdRab(1) += dEtuv(1) * Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dQdRab(2) += dEtuv(2) * Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));

                        }
                    }
                }

            }
        }
    }

    return dQdRab;

}

rowvec Integrator::electronRepulsionIntegral_Rcd_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        int iC, int jC, int kC, int iD, int jD, int kD)
{

    rowvec dQdRcd = zeros<rowvec>(3);

    int tMax = iA + iB + 1;
    int uMax = jA + jB + 1;
    int vMax = kA + kB + 1;
    int kMax = iC + iD + 1;
    int lMax = jC + jD + 1;
    int mMax = kC + kD + 1;

    for(int t = 0; t < tMax; t++){
        for(int u = 0; u < uMax; u++){
            for(int v = 0; v < vMax; v++){

                double Etuv= m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v);

                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            rowvec3 dEklm = zeros<rowvec>(3);;
                            dEklm(0) = m_dEcd(0)(iC, iD, t) * m_Ecd(1)(jC, jD, u) * m_Ecd(2)(kC, kD, v);
                            dEklm(1) = m_Ecd(0)(iC, iD, t) * m_dEcd(1)(jC, jD, u) * m_Ecd(2)(kC, kD, v);
                            dEklm(2) = m_Ecd(0)(iC, iD, t) * m_Ecd(1)(jC, jD, u) * m_dEcd(2)(kC, kD, v);

                            dQdRcd(0) += dEklm(0) * Etuv * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;
                            dQdRcd(1) += dEklm(1) * Etuv * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;
                            dQdRcd(2) += dEklm(2) * Etuv * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;

                        }
                    }
                }

            }
        }
    }

    return dQdRcd;

}










