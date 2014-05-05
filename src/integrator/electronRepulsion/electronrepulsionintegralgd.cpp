#include "electronrepulsionintegralgd.h"

using namespace hf;
ElectronRepulsionIntegralGD::ElectronRepulsionIntegralGD(const int highestOrder,
                                                         const field<cube> *Eab,
                                                         const field<cube> *Ecd,
                                                         const field<cube> *dEab_dQab,
                                                         const field<cube> *dEcd_dQcd,
                                                         const PrimitiveGTO *primitiveA,
                                                         const PrimitiveGTO *primitiveB,
                                                         const PrimitiveGTO *primitiveC,
                                                         const PrimitiveGTO *primitiveD):
    m_R(new HermiteIntegrals(highestOrder)),
    m_Eab(Eab),
    m_Ecd(Ecd),
    m_dEab_dQab(dEab_dQab),
    m_dEcd_dQcd(dEcd_dQcd),
    m_primitiveA(primitiveA),
    m_primitiveB(primitiveB),
    m_primitiveC(primitiveC),
    m_primitiveD(primitiveD)
{
}

void ElectronRepulsionIntegralGD::updateHermiteIntegrals()
{

    const rowvec &A = m_primitiveA->center();
    const rowvec &B = m_primitiveB->center();
    const rowvec &C = m_primitiveC->center();
    const rowvec &D = m_primitiveD->center();

    const double &a  = m_primitiveA->exponent();
    const double &b  = m_primitiveB->exponent();
    const double &c  = m_primitiveC->exponent();
    const double &d  = m_primitiveD->exponent();

    double p = a + b;
    double q = c + d;
    double alpha = p*q/(p+q);
    rowvec PQ = (a*A + b*B)/p - (c*C + d*D)/q;

    m_tMax = m_primitiveA->xPower() + m_primitiveB->xPower() + 1;
    m_uMax = m_primitiveA->yPower() + m_primitiveB->yPower() + 1;
    m_vMax = m_primitiveA->zPower() + m_primitiveB->zPower() + 1;
    m_kMax = m_primitiveC->xPower() + m_primitiveD->xPower() + 1;
    m_lMax = m_primitiveC->yPower() + m_primitiveD->yPower() + 1;
    m_mMax = m_primitiveC->zPower() + m_primitiveD->zPower() + 1;
    m_R->updateR(PQ, alpha,
                 m_tMax + m_kMax - 1, m_uMax + m_lMax - 1, m_vMax + m_mMax - 1);

}

rowvec ElectronRepulsionIntegralGD::QabDerivative()
{
    rowvec dGdQab = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();
    int iC = m_primitiveC->xPower();
    int jC = m_primitiveC->yPower();
    int kC = m_primitiveC->zPower();
    int iD = m_primitiveD->xPower();
    int jD = m_primitiveD->yPower();
    int kD = m_primitiveD->zPower();

    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){

                rowvec3 dEtuv = zeros<rowvec>(3);
                dEtuv(0) = m_dEab_dQab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v);
                dEtuv(1) = m_Eab->at(0)(iA, iB, t) * m_dEab_dQab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v);
                dEtuv(2) = m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_dEab_dQab->at(2)(kA, kB, v);

                for(int k = 0; k < m_kMax; k++){
                    for(int l = 0; l < m_lMax; l++){
                        for(int m = 0; m < m_mMax; m++){

                            double Eklm= m_Ecd->at(0)(iC, iD, k) * m_Ecd->at(1)(jC, jD, l) * m_Ecd->at(2)(kC, kD, m);

                            dGdQab(0) += dEtuv(0) * Eklm * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dGdQab(1) += dEtuv(1) * Eklm * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dGdQab(2) += dEtuv(2) * Eklm * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));

                        }
                    }
                }

            }
        }
    }

    return dGdQab * m_primitiveA->weight() * m_primitiveB->weight()
                  * m_primitiveC->weight() * m_primitiveD->weight();

}

rowvec ElectronRepulsionIntegralGD::QcdDerivative()
{
    rowvec dGdQcd = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();
    int iC = m_primitiveC->xPower();
    int jC = m_primitiveC->yPower();
    int kC = m_primitiveC->zPower();
    int iD = m_primitiveD->xPower();
    int jD = m_primitiveD->yPower();
    int kD = m_primitiveD->zPower();

    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){

                double Etuv= m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v);

                for(int k = 0; k < m_kMax; k++){
                    for(int l = 0; l < m_lMax; l++){
                        for(int m = 0; m < m_mMax; m++){

                            rowvec3 dEklm = zeros<rowvec>(3);
                            dEklm(0) = m_dEcd_dQcd->at(0)(iC, iD, k) * m_Ecd->at(1)(jC, jD, l) * m_Ecd->at(2)(kC, kD, m);
                            dEklm(1) = m_Ecd->at(0)(iC, iD, k) * m_dEcd_dQcd->at(1)(jC, jD, l) * m_Ecd->at(2)(kC, kD, m);
                            dEklm(2) = m_Ecd->at(0)(iC, iD, k) * m_Ecd->at(1)(jC, jD, l) * m_dEcd_dQcd->at(2)(kC, kD, m);

                            dGdQcd(0) += dEklm(0) * Etuv * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;
                            dGdQcd(1) += dEklm(1) * Etuv * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;
                            dGdQcd(2) += dEklm(2) * Etuv * m_R->R(0,t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));;

                        }
                    }
                }

            }
        }
    }

    return dGdQcd* m_primitiveA->weight() * m_primitiveB->weight()
            * m_primitiveC->weight() * m_primitiveD->weight();

}




rowvec ElectronRepulsionIntegralGD::PabDerivative()
{
    rowvec dGdPab = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();
    int iC = m_primitiveC->xPower();
    int jC = m_primitiveC->yPower();
    int kC = m_primitiveC->zPower();
    int iD = m_primitiveD->xPower();
    int jD = m_primitiveD->yPower();
    int kD = m_primitiveD->zPower();

    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){

                double Etuv= m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v);

                for(int k = 0; k < m_kMax; k++){
                    for(int l = 0; l < m_lMax; l++){
                        for(int m = 0; m < m_mMax; m++){

                            double Eklm= m_Ecd->at(0)(iC, iD, k) * m_Ecd->at(1)(jC, jD, l) * m_Ecd->at(2)(kC, kD, m);

                            dGdPab(0) += Etuv*Eklm  * m_R->R(0,t+1+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dGdPab(1) += Etuv*Eklm  * m_R->R(0,t+k,u+1+l,v+m) * (1 - 2* ((k+l+m)%2));
                            dGdPab(2) += Etuv*Eklm  * m_R->R(0,t+k,u+l,v+1+m) * (1 - 2* ((k+l+m)%2));

                        }
                    }
                }

            }
        }
    }

    return dGdPab* m_primitiveA->weight() * m_primitiveB->weight()
            * m_primitiveC->weight() * m_primitiveD->weight();

}


rowvec ElectronRepulsionIntegralGD::PcdDerivative()
{
    rowvec dGdPcd = zeros<rowvec>(3);

    int iA = m_primitiveA->xPower();
    int jA = m_primitiveA->yPower();
    int kA = m_primitiveA->zPower();
    int iB = m_primitiveB->xPower();
    int jB = m_primitiveB->yPower();
    int kB = m_primitiveB->zPower();
    int iC = m_primitiveC->xPower();
    int jC = m_primitiveC->yPower();
    int kC = m_primitiveC->zPower();
    int iD = m_primitiveD->xPower();
    int jD = m_primitiveD->yPower();
    int kD = m_primitiveD->zPower();

    for(int t = 0; t < m_tMax; t++){
        for(int u = 0; u < m_uMax; u++){
            for(int v = 0; v < m_vMax; v++){

                double Etuv= m_Eab->at(0)(iA, iB, t) * m_Eab->at(1)(jA, jB, u) * m_Eab->at(2)(kA, kB, v);

                for(int k = 0; k < m_kMax; k++){
                    for(int l = 0; l < m_lMax; l++){
                        for(int m = 0; m < m_mMax; m++){

                            double Eklm= m_Ecd->at(0)(iC, iD, k) * m_Ecd->at(1)(jC, jD, l) * m_Ecd->at(2)(kC, kD, m);

                            dGdPcd(0) += Etuv*Eklm  * m_R->R(0, t+k+1,u+l,v+m) * (1 - 2* ((k+1+l+m)%2));
                            dGdPcd(1) += Etuv*Eklm  * m_R->R(0, t+k,u+l+1,v+m) * (1 - 2* ((k+l+1+m)%2));
                            dGdPcd(2) += Etuv*Eklm  * m_R->R(0, t+k,u+l,v+m+1) * (1 - 2* ((k+l+m+1)%2));

                        }
                    }
                }

            }
        }
    }

    return dGdPcd* m_primitiveA->weight() * m_primitiveB->weight()
            * m_primitiveC->weight() * m_primitiveD->weight();

}

