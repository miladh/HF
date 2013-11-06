#include "hermitecoefficients.h"

HermiteCoefficients::HermiteCoefficients()
{
}

bool HermiteCoefficients::interiorPoint(int iA, int iB, int t)
{
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}


void HermiteCoefficients::setupE(const double &a, const rowvec3 &A, const int &La,
                                 const double &b, const rowvec3 &B, const int &Lb,
                                 field<cube> &E)
{
    int iAmax = La + 3;
    int iBmax = Lb + 3;
    int tmax  = iAmax + iBmax - 1;

    for(int cor = 0; cor < 3; cor++){
        E(cor) = zeros(iAmax, iBmax, tmax);
    }

    double p = a + b;
    rowvec P  = (a*A + b*B)/p;
    rowvec AB = A - B;
    rowvec PA = P - A;
    rowvec PB = P - B;
    rowvec Kab = exp(-(a * b / p)*AB%AB);

    for(uint cor = 0; cor < E.n_elem; cor++){
        E[cor](0,0,0) = Kab(cor);
    }

    for(uint cor=0; cor < E.n_elem; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
        for(int iB = 1; iB < iBmax; iB++){
            for(int t = 0; t < tmax; t++){

                int iA = 0;
                int iBp = iB - 1;
                int tp = t - 1;
                int tn = t + 1;

                double E_iA_iBp_tp = 0.0;
                if(interiorPoint(iA, iBp, tp)){
                    E_iA_iBp_tp = E[cor](iA, iBp, tp);
                }

                double E_iA_iBp_t = 0;
                if(interiorPoint(iA, iBp, t)) {
                    E_iA_iBp_t = E[cor](iA, iBp, t);
                }

                double E_iA_iBp_tn = 0;
                if(interiorPoint(iA, iBp, tn)) {
                    E_iA_iBp_tn = E[cor](iA, iBp, tn);
                }

                E[cor](iA,iB,t) = 1.0 / (2*p) * E_iA_iBp_tp + PB(cor) * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
            }
        }



        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
        for(int iA = 1; iA < iAmax; iA++) {
            for(int iB = 0; iB < iBmax; iB++) {
                for(int t = 0; t < tmax; t++) {

                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;

                    double E_iAp_iB_tp = 0;
                    if(interiorPoint(iAp, iB, tp)) {
                        E_iAp_iB_tp = E[cor](iAp, iB, tp);
                    }

                    double E_iAp_iB_t = 0;
                    if(interiorPoint(iAp, iB, t)) {
                        E_iAp_iB_t = E[cor](iAp, iB, t);
                    }

                    double E_iAp_iB_tn = 0;
                    if(interiorPoint(iAp, iB, tn)) {
                        E_iAp_iB_tn = E[cor](iAp, iB, tn);
                    }

                    E[cor](iA,iB,t) = 1.0 / (2*p) * E_iAp_iB_tp + PA(cor) * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }

    }//End of cor=(x,y,z) loop
}
