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


    Eab = new HermiteCoefficients(maxAngularMomentum);
    Ecd = new HermiteCoefficients(maxAngularMomentum);

    m_overlap = new OverlapIntegral(Eab->coefficients(), &m_primitiveA, &m_primitiveB);
    m_kinetic = new KineticIntegral(m_overlap, &m_primitiveA, &m_primitiveB);
    m_nuclearAttraction = new NuclearAttractionIntegral(2 * maxAngularMomentum + 1,
                                                        Eab->coefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_corePositionC);

    m_electronRepulsion = new ElectronRepulsionIntegral(4 * maxAngularMomentum + 1,
                                                        Eab->coefficients(),
                                                        Ecd->coefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_primitiveC,
                                                        &m_primitiveD);



    m_overlapGD = new OverlapIntegralGD(m_overlap, Eab->QDerivativeCoefficients());
    m_kineticGD = new KineticIntegralGD(m_kinetic, m_overlapGD);
    m_nuclearAttractionGD = new NuclearAttractionIntegralGD(2 * maxAngularMomentum + 1,
                                                        Eab->coefficients(),
                                                        Eab->QDerivativeCoefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_corePositionC);




//    cout << "--------------"<< endl;
//    cout << "Integrator: "  << endl;
//    cout << "primA adr   "  << &m_primitiveA.center()<< endl;
//    cout << "primB adr   "  << &m_primitiveA<< endl;
//    sleep(4);
//exit(0);
}

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
            Eab->updateE(m_primitiveA, m_primitiveB);
            m_hermiteCoefficients.setupE(m_primitiveA, m_primitiveB, m_Eab);
        }else{
            Eab->updateE(m_primitiveA, m_primitiveB,false);
            m_hermiteCoefficients.setupE(m_primitiveA, m_primitiveB, m_Eab,false);
        }

    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setupE(m_primitiveC, m_primitiveD, m_Ecd, false);
        Ecd->updateE(m_primitiveC, m_primitiveD, false);

    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }
}

void Integrator::updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral,bool kin)
{

    if(oneParticleIntegral){
        if(kin){
        m_hermiteCoefficients.setup_dEdR(m_primitiveA, m_primitiveB, m_Eab, m_dEab);
        Eab->updatedE_dQ(m_primitiveA, m_primitiveB);
        }else{
            m_hermiteCoefficients.setup_dEdR(m_primitiveA, m_primitiveB, m_Eab, m_dEab, false);
            Eab->updatedE_dQ(m_primitiveA, m_primitiveB,false);

        }

    }else if(twoParticleIntegral){
        m_hermiteCoefficients.setup_dEdR(m_primitiveC, m_primitiveD, m_Ecd,m_dEcd, false);
        Ecd->updatedE_dQ(m_primitiveA, m_primitiveB, false);

    }else{
        cerr << "Hermite coefficients not updated!" << endl;
    }

}

void Integrator::updateHermiteIntegrals()
{
    m_nuclearAttractionGD->updateHermiteIntegrals();
}
/********************************************************************************************
 *
 *                                  Molecular Gaussian Integrals
 *
 * ******************************************************************************************/
double Integrator::overlapIntegral()
{
    return m_overlap->evaluate();
}

double Integrator::kineticIntegral()
{
    return m_kinetic->evaluate();
}

double Integrator::nuclearAttractionIntegral()
{
    return m_nuclearAttraction->evaluate();
}

double Integrator::electronRepulsionIntegral()
{
    return m_electronRepulsion->evaluate();
}


/********************************************************************************************
 *
 *                  Molecular Gaussian Integral Geometrical Derivatives (GD)
 *
 * ******************************************************************************************/

rowvec Integrator::QDerivativeOverlapIntegral()
{
    return m_overlapGD->evaluate() * m_primitiveA.weight() * m_primitiveB.weight();
}


rowvec Integrator::QDerivativeKineticIntegral() {

    return m_kineticGD->evaluate() * m_primitiveA.weight() * m_primitiveB.weight();
}

rowvec Integrator::QDerivativeNuclearAttractionIntegral()
{

    return m_nuclearAttractionGD->QDerivative();
}

rowvec Integrator::PDerivativeNuclearAttractionIntegral()
{

    return m_nuclearAttractionGD->PDerivative();
}

rowvec Integrator::CDerivativeNuclearAttractionIntegral()
{
    return m_nuclearAttractionGD->CDerivative();
}

/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------*/
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










