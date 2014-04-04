#include "hermitecoefficients.h"
using namespace hf;
HermiteCoefficients::HermiteCoefficients()
{
}

HermiteCoefficients::HermiteCoefficients(const int maxAngularMomentum)
{
    m_E.set_size(3);
    m_dE_dQ.set_size(3);

    int iAmax = maxAngularMomentum + 3;
    int iBmax = maxAngularMomentum + 3;
    int tmax  = iAmax + iBmax - 1;

    for(int cor = 0; cor < 3; cor++){
        m_E(cor)  =zeros(iAmax, iBmax, tmax);
    }

    for(int cor = 0; cor < 3; cor++){
        m_dE_dQ(cor) = zeros(iAmax , iBmax, tmax);
    }
}

bool HermiteCoefficients::interiorPoint(int iA, int iB, int t)
{
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}
const field<cube>* HermiteCoefficients::coefficients()const
{
    return &m_E;
}


const field<cube>* HermiteCoefficients::QDerivativeCoefficients()const
{
    return &m_dE_dQ;
}

void HermiteCoefficients::updateE(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB,
                                  bool kin)
{
    const rowvec R = primitiveA.center() - primitiveB.center();
    rowvec iAmax = primitiveA.powers() + 1;
    rowvec iBmax = primitiveB.powers() + 1;

    if(kin){
        iAmax += 2;
        iBmax += 2;
    }
    rowvec tmax  = iAmax + iBmax - 1;

    double a = primitiveA.exponent();
    double b = primitiveB.exponent();
    double p = a + b;
    double factor = -(a * b / p);

    for(uint cor = 0; cor < m_E.n_elem; cor++){
        m_E(cor)(0,0,0) = std::exp(factor*R(cor)*R(cor));
    }

    for(uint cor=0; cor < m_E.n_elem; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPB * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
        for(int iB = 1; iB < iBmax(cor); iB++){
            for(int t = 0; t < tmax(cor); t++){

                int iA = 0;
                int iBp = iB - 1;
                int tp = t - 1;
                int tn = t + 1;

                double E_iA_iBp_tp = 0.0;
                if(interiorPoint(iA, iBp, tp)){
                    E_iA_iBp_tp = m_E(cor)(iA, iBp, tp);
                }

                double E_iA_iBp_t = 0;
                if(interiorPoint(iA, iBp, t)) {
                    E_iA_iBp_t = m_E(cor)(iA, iBp, t);
                }

                double E_iA_iBp_tn = 0;
                if(interiorPoint(iA, iBp, tn)) {
                    E_iA_iBp_tn = m_E(cor)(iA, iBp, tn);
                }

                m_E(cor)(iA,iB,t) = 1.0 / (2*p) * E_iA_iBp_tp + a / p* R(cor) * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
            }
        }



        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
        for(int iA = 1; iA < iAmax(cor); iA++) {
            for(int iB = 0; iB < iBmax(cor); iB++) {
                for(int t = 0; t < tmax(cor); t++) {

                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;

                    double E_iAp_iB_tp = 0;
                    if(interiorPoint(iAp, iB, tp)) {
                        E_iAp_iB_tp = m_E(cor)(iAp, iB, tp);
                    }

                    double E_iAp_iB_t = 0;
                    if(interiorPoint(iAp, iB, t)) {
                        E_iAp_iB_t = m_E(cor)(iAp, iB, t);
                    }

                    double E_iAp_iB_tn = 0;
                    if(interiorPoint(iAp, iB, tn)) {
                        E_iAp_iB_tn = m_E(cor)(iAp, iB, tn);
                    }

                    m_E(cor)(iA,iB,t) = 1.0 / (2*p) * E_iAp_iB_tp - b / p* R(cor) * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }

    }//End of cor=(x,y,z) loop
}





void HermiteCoefficients::updatedE_dQ(const PrimitiveGTO &primitiveA, const PrimitiveGTO &primitiveB, bool kin)
{

    const rowvec R = primitiveA.center() - primitiveB.center();
    rowvec iAmax = primitiveA.powers() + 1;
    rowvec iBmax = primitiveB.powers() + 1;

    if(kin){
        iAmax += 2;
        iBmax += 2;
    }
    rowvec tmax  = iAmax + iBmax - 1;

    double a = primitiveA.exponent();
    double b = primitiveB.exponent();
    double p = a + b;
    double factor = -2 * a * b / p;

    for(uint cor = 0; cor < m_dE_dQ.n_elem; cor++){
        m_dE_dQ(cor)(0,0,0) = factor * R(cor) * m_E(cor)(0,0,0);
    }

    for(uint cor=0; cor < m_dE_dQ.n_elem; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // dE(t,i,j) = 1 / (2*p) * dE(t-1,i,j-1) + a/p ( Rcor * dE(t,i,j-1) + E(t,i,j-1) ) + (t + 1)*dE(t+1,i,j-1)
        for(int iB = 1; iB < iBmax(cor); iB++){
            for(int t = 0; t < tmax(cor); t++){

                int iA = 0;
                int iBp = iB - 1;
                int tp = t - 1;
                int tn = t + 1;

                double dE_iA_iBp_tp = 0.0;
                if(interiorPoint(iA, iBp, tp)){
                    dE_iA_iBp_tp = m_dE_dQ(cor)(iA, iBp, tp);
                }

                double dE_iA_iBp_t = 0;
                double E_iA_iBp_t  = 0;
                if(interiorPoint(iA, iBp, t)) {
                    dE_iA_iBp_t = m_dE_dQ(cor)(iA, iBp, t);
                    E_iA_iBp_t = m_E(cor)(iA, iBp, t);
                }

                double dE_iA_iBp_tn = 0;
                if(interiorPoint(iA, iBp, tn)) {
                    dE_iA_iBp_tn = m_dE_dQ(cor)(iA, iBp, tn);
                }

                m_dE_dQ(cor)(iA,iB,t) = 1.0 / (2*p) * dE_iA_iBp_tp + a / p * (R(cor)* dE_iA_iBp_t + E_iA_iBp_t)
                        + (t + 1)*dE_iA_iBp_tn;
            }
        }



        // p = previous, n = next
        // dE(t,i,j) = 1 / (2*p) * dE(t-1,i-1,j) - b/p ( Rcor * dE(t,i-1,j) + E(t,i-1,j) ) + (t + 1)*dE(t+1,i-1,j)
        for(int iA = 1; iA < iAmax(cor); iA++) {
            for(int iB = 0; iB < iBmax(cor); iB++) {
                for(int t = 0; t < tmax(cor); t++) {

                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;

                    double dE_iAp_iB_tp = 0;
                    if(interiorPoint(iAp, iB, tp)) {
                        dE_iAp_iB_tp = m_dE_dQ(cor)(iAp, iB, tp);
                    }

                    double dE_iAp_iB_t = 0;
                    double E_iAp_iB_t = 0;
                    if(interiorPoint(iAp, iB, t)) {
                        dE_iAp_iB_t = m_dE_dQ(cor)(iAp, iB, t);
                        E_iAp_iB_t = m_E(cor)(iAp, iB, t);
                    }

                    double dE_iAp_iB_tn = 0;
                    if(interiorPoint(iAp, iB, tn)) {
                        dE_iAp_iB_tn = m_dE_dQ(cor)(iAp, iB, tn);
                    }

                    m_dE_dQ(cor)(iA,iB,t) = 1.0 / (2*p) * dE_iAp_iB_tp - b / p * (R(cor) *dE_iAp_iB_t  + E_iAp_iB_t )
                            +  (t + 1)*dE_iAp_iB_tn;
                }
            }
        }

    }//End of cor=(x,y,z) loop
}




































