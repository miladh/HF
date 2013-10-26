#include "integrator.h"

Integrator::Integrator():
    m_corePositionA(zeros<rowvec>(3)),
    m_corePositionB(zeros<rowvec>(3))
{
    setMaxAngularMomentum(0);
}


void Integrator::setCorePositions(const mat corePositions){

    m_corePositions = corePositions;
}

mat Integrator::corePositions() const
{
    return m_corePositions;
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

uint Integrator::maxAngularMomentum() const
{
    return m_maxAngularMomentum;
}

void Integrator::setMaxAngularMomentum(const uint &maxAngularMomentum)
{
    m_maxAngularMomentum = maxAngularMomentum;
}

void Integrator::addPrimitives(PrimitiveGTO *primitive)
{
    m_primitives.push_back(primitive);

}

bool Integrator::interiorPoint(int iA, int iB, int t)
{
    if(t < 0 || t > (iA + iB) || iA < 0 || iB < 0) {
        return false;
    } else {
        return true;
    }
}


void Integrator::setupR(const rowvec &C){

//    uint iAmax = m_maxAngularMomentum + 3;
//    uint iBmax = m_maxAngularMomentum + 3;
//    uint tmax  = 2*(iAmax +iBmax);

//    for(uint n = 0; n < 6*tmax; n++){
//        m_R.push_back(zeros(tmax, tmax, tmax));
//    }


//    uint nMax = 3 * tmax;
//    const rowvec &A = m_corePositionA;
//    const rowvec &B = m_corePositionB;

//    double a = m_primitives[0]->exponent();
//    double b = m_primitives[1]->exponent();

//    double p = a + b;
//    rowvec P  = (a*A + b*B)/p;
//    rowvec PC = P - C;
//    Boys boys(p*dot(PC,PC),nMax);
//    rowvec Fn = boys.getBoysFunctions();


//    for(uint n = 0; n < nMax; n++){
//        m_R.at(n)(0,0,0) = pow(-2*p,n)*Fn[n];
//    }

//    for(uint tuvSum = 1; tuvSum <= nMax; tuvSum++){
//        for(uint n = 0; n <= tuvSum; n++){

//            for(uint t = 0; t < tuvSum; t++){
//                for(uint u = 0; u < tuvSum; u++){
//                    for(uint v = 0; v < tuvSum; v++){
//                        if(t+u+v == tuvSum){

//                            m_R.at(n)(t+1,u,v) = t*m_R.at(n+1)(t-1,u,v) + PC(0) * m_R.at(n+1)(t,u,v);
//                            m_R.at(n)(t,u+1,v) = u*m_R.at(n+1)(t,u-1,v) + PC(1) * m_R.at(n+1)(t,u,v);
//                            m_R.at(n)(t,u,v+1) = v*m_R.at(n+1)(t,u,v-1) + PC(2) * m_R.at(n+1)(t,u,v);
//                        }
//                    }
//                }
//            }
//        }
//    }

}


void Integrator::setupE()
{
    int iAmax = m_maxAngularMomentum + 3;
    int iBmax = m_maxAngularMomentum + 3;
    int tmax  = 2*(iAmax +iBmax);

    for(int cor = 0; cor < 3; cor++){
        m_E[cor] = zeros(iAmax, iBmax, tmax);
    }

    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;

    double a = m_primitives[0]->exponent();
    double b = m_primitives[1]->exponent();

    double p = a + b;
    double mu = a * b / p;

    rowvec P  = (a*A + b*B)/p;
    rowvec AB = A - B;
    rowvec PA = P - A;
    rowvec PB = P - B;
    rowvec Kab = exp(-mu*AB%AB);


    for(int cor = 0; cor < 3; cor++){
        m_E[cor](0,0,0) = Kab(cor);
    }

    for(int cor=0; cor < 3; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
        for(int iB = 1; iB < iBmax; iB++){
            for(int t = 0; t < tmax-1; t++){

                int iA = 0;
                int iBp = iB - 1;
                int tp = t - 1;
                int tn = t + 1;

                double E_iA_iBp_tp = 0.0;
                if(interiorPoint(iA, iBp, tp)){
                    E_iA_iBp_tp = m_E[cor](iA, iBp, tp);
                }

                double E_iA_iBp_t = 0;
                if(interiorPoint(iA, iBp, t)) {
                    E_iA_iBp_t = m_E[cor](iA, iBp, t);
                }

                double E_iA_iBp_tn = 0;
                if(interiorPoint(iA, iBp, tn)) {
                    E_iA_iBp_tn = m_E[cor](iA, iBp, tn);
                }

                m_E[cor](iA,iB,t) = 1 / (2*p) * E_iA_iBp_tp + PB(cor) * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
            }
        }



        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i-1,j) + XPA * E(t,i-1,j) + (t + 1)*E(t+1,i-1,j)
        for(int iA = 1; iA < iAmax; iA++) {
            for(int iB = 0; iB < iBmax; iB++) {
                for(int t = 0; t < tmax - 1; t++) {

                    int iAp = iA - 1;
                    int tp = t - 1;
                    int tn = t + 1;

                    double E_iAp_iB_tp = 0;
                    if(interiorPoint(iAp, iB, tp)) {
                        E_iAp_iB_tp = m_E[cor](iAp, iB, tp);
                    }

                    double E_iAp_iB_t = 0;
                    if(interiorPoint(iAp, iB, t)) {
                        E_iAp_iB_t = m_E[cor](iAp, iB, t);
                    }

                    double E_iAp_iB_tn = 0;
                    if(interiorPoint(iAp, iB, tn)) {
                        E_iAp_iB_tn = m_E[cor](iAp, iB, tn);
                    }

                    m_E[cor](iA,iB,t) = 1 / (2*p) * E_iAp_iB_tp + PA(cor) * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
                }
            }
        }

    }//End of cor=(x,y,z) loop


    //    cout <<"p " << p << endl;
    //    cout <<"mu "<< mu << endl;
    //    cout <<"P"<< P << endl;
    //    cout <<"ab"<< AB <<endl;
    //    cout <<"pa"<< PA <<endl;
    //    cout <<"pb"<< PB <<endl;
    //    cout <<"kab"<< Kab <<endl;


    //    cout << m_E[0] << endl;
    //    cout << "-----------------------------------" <<endl;
    //    cout << m_E[1] << endl;
    //    cout << "-----------------------------------" <<endl;
    //    cout << m_E[2] << endl;

}


double Integrator::overlapIntegral(int cor, int iA, int iB)
{
    double a = m_primitives[0]->exponent();
    double b = m_primitives[1]->exponent();
    double p = a + b;
    return m_E[cor](iA,iB,0) * sqrt(M_PI / p);
}


double Integrator::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {

    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}


double Integrator::kineticIntegral(int dim, int iA, int iB) {
    double b = m_primitives[1]->exponent();

    double S_iA_iBnn = overlapIntegral(dim, iA, iB + 2);
    double S_iA_iB = overlapIntegral(dim, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= overlapIntegral(dim, iA, iB - 2);
    } else {
        S_iA_iBpp = 0;
    }
    return 4 * b * b * S_iA_iBnn - 2*b * (2*iB + 1) * S_iA_iB + iB * (iB + 1) * S_iA_iBpp;
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













