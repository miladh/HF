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
    //        setMaxAngularMomentum(0);
}

void Integrator::setMaxAngularMomentum(const uint &maxAngularMomentum)
{
    m_maxAngularMomentum = maxAngularMomentum;


    m_Eab.set_size(3);  m_Ecd.set_size(3);
    m_dEab.set_size(3); m_dEcd.set_size(3);
    int iAmax = m_maxAngularMomentum + 3;
    int iBmax = m_maxAngularMomentum + 3;
    int tmax  = iAmax + iBmax - 1;
    for(int cor = 0; cor < 3; cor++){
        m_Eab(cor)  = zeros(iAmax, iBmax, tmax);
        m_Ecd(cor)  = zeros(iAmax, iBmax, tmax);
        m_dEab(cor) = zeros(iAmax, iBmax, tmax);
        m_dEcd(cor) = zeros(iAmax, iBmax, tmax);
    }

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


void Integrator::updateHermiteCoefficients(bool oneParticleIntegral, bool twoParticleIntegral)
{

    if(oneParticleIntegral){
        m_hermiteCoefficients.setupE(m_exponentA, m_exponentB,
                                     m_corePositionA - m_corePositionB,
                                     m_Eab);

    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setupE(m_exponentC, m_exponentD,
                                     m_corePositionC - m_corePositionD,
                                     m_Ecd);
    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }
}

void Integrator::updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral)
{

    if(oneParticleIntegral){
        m_hermiteCoefficients.setup_dEdR(m_exponentA, m_exponentB,
                                         m_corePositionA - m_corePositionB,
                                         m_Eab, m_dEab);


    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setup_dEdR(m_exponentC, m_exponentD,
                                         m_corePositionC - m_corePositionD,
                                         m_Ecd,m_dEcd);
    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }

}


/*---------------------------------------------------------------------------------------------------*/
double Integrator::overlapIntegral(int cor, int iA, int iB)
{
    double p = m_exponentA + m_exponentB;
    return m_Eab[cor](iA,iB,0) * sqrt(M_PI / p);
}


double Integrator::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB)
{

    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}


double Integrator::overlapIntegral_derivative(int cor, int iA, int iB)
{
    double p =  m_exponentA +  m_exponentB;
    return m_dEab(cor)(iA,iB,0) * sqrt(M_PI / p);
}


rowvec Integrator::overlapIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB)
{
    rowvec dS = zeros<rowvec>(3);

    dS(0) = overlapIntegral_derivative(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
    dS(1) = overlapIntegral_derivative(1, jA, jB) * overlapIntegral(0, iA, iB) * overlapIntegral(2, kA, kB);
    dS(2) = overlapIntegral_derivative(2, kA, kB) * overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB);

    return dS;
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



double Integrator::kineticIntegral_derivative(int cor, int iA, int iB) {
    double b = m_exponentB;

    double dS_iA_iBnn = overlapIntegral_derivative(cor, iA, iB + 2);
    double dS_iA_iB = overlapIntegral_derivative(cor, iA, iB);
    double dS_iA_iBpp;

    if(iB - 2 >= 0) {
        dS_iA_iBpp= overlapIntegral_derivative(cor, iA, iB - 2);
    } else {
        dS_iA_iBpp = 0;
    }
    return 4 * b * b * dS_iA_iBnn - 2*b * (2*iB + 1) * dS_iA_iB + iB * (iB - 1) * dS_iA_iBpp;
}


rowvec Integrator::kineticIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB) {
    rowvec dT = zeros<rowvec>(3);

    double T_iA_iB = kineticIntegral(0, iA, iB);
    double T_jA_jB = kineticIntegral(1, jA, jB);
    double T_kA_kB = kineticIntegral(2, kA, kB);

    double dT_iA_iB = kineticIntegral_derivative(0, iA, iB);
    double dT_jA_jB = kineticIntegral_derivative(1, jA, jB);
    double dT_kA_kB = kineticIntegral_derivative(2, kA, kB);

    double S_iA_iB = overlapIntegral(0, iA, iB);
    double S_jA_jB = overlapIntegral(1, jA, jB);
    double S_kA_kB = overlapIntegral(2, kA, kB);

    double dS_iA_iB = overlapIntegral_derivative(0, iA, iB);
    double dS_jA_jB = overlapIntegral_derivative(1, jA, jB);
    double dS_kA_kB = overlapIntegral_derivative(2, kA, kB);

    dT(0) = (dT_iA_iB * S_jA_jB * S_kA_kB) + (dS_iA_iB * T_jA_jB * S_kA_kB) + (dS_iA_iB * S_jA_jB  * T_kA_kB);
    dT(1) = (T_iA_iB * dS_jA_jB * S_kA_kB) + (S_iA_iB * dT_jA_jB * S_kA_kB) + (S_iA_iB * dS_jA_jB  * T_kA_kB);
    dT(2) = (T_iA_iB * S_jA_jB * dS_kA_kB) + (S_iA_iB * T_jA_jB * dS_kA_kB) + (S_iA_iB * S_jA_jB  * dT_kA_kB);

    dT *= -0.5;
    return dT;
}


/*---------------------------------------------------------------------------------------------------*/
double Integrator::nuclearAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB)
{
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;

    const double &a  = m_exponentA;
    const double &b  = m_exponentB;

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
                result += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v);
            }
        }
    }

    return 2 * result * M_PI / p;
}



rowvec Integrator::nuclearAttractionIntegral_R_derivative(int iA, int jA, int kA, int iB, int jB, int kB)
{
    rowvec dVdR = zeros<rowvec>(3);

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
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

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
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

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;


    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
                dVdC(0) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t+1,u,v);
                dVdC(1) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u+1,v);
                dVdC(2) += m_Eab(0)(iA, iB, t) * m_Eab(1)(jA, jB, u) * m_Eab(2)(kA, kB, v) * m_Ren(0)(t,u,v+1);
            }
        }
    }

    return dVdC;
}

rowvec Integrator::nuclearAttractionIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        bool differentiateWrtA, bool differentiateWrtB,
                                                        bool differentiateWrtC)
{
    rowvec dVab = zeros<rowvec>(3);
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;

    const double &a  = m_exponentA;
    const double &b  = m_exponentB;

    double p = a + b;
    rowvec P  = (a*A + b*B)/p;
    rowvec PC = P - C;

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

double Integrator::electronRepulsionIntegral(int iA, int jA, int kA, int iB, int jB, int kB,
                                             int iC, int jC, int kC, int iD, int jD, int kD)
{
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;
    const rowvec &D = m_corePositionD;

    const double &a  = m_exponentA;
    const double &b  = m_exponentB;
    const double &c  = m_exponentC;
    const double &d  = m_exponentD;

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
                            result += Etuv*Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                        }
                    }
                }

            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q));

    return result;

}
rowvec Integrator::electronRepulsionIntegral_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                        int iC, int jC, int kC, int iD, int jD, int kD,
                                                        bool differentiateWrtA, bool differentiateWrtB,
                                                        bool differentiateWrtC, bool differentiateWrtD)
{
    rowvec dQabcd = zeros<rowvec>(3);
    rowvec dQdPab, dQdPcd;
    rowvec dQdRab, dQdRcd;
    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;
    const rowvec &D = m_corePositionD;

    const double &a  = m_exponentA;
    const double &b  = m_exponentB;
    const double &c  = m_exponentC;
    const double &d  = m_exponentD;

    double p = a + b;
    double q = c + d;
    rowvec P  = (a*A + b*B)/p;
    rowvec Q  = (c*C + d*D)/q;

    double alpha = p*q/(p+q);
    rowvec PQ = P - Q;

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










