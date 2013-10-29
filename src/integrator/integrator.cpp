#include "integrator.h"

Integrator::Integrator():
    m_corePositionA(zeros<rowvec>(3)),
    m_corePositionB(zeros<rowvec>(3))
{
    setMaxAngularMomentum(0);
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


void Integrator::setupR(const rowvec &PQ, const double &alpha){

    uint tMax  = 2 * m_maxAngularMomentum + 1;
    uint uMax  = 2 * m_maxAngularMomentum + 1;
    uint vMax  = 2 * m_maxAngularMomentum + 1;
    uint nMax  = 2 * m_maxAngularMomentum + 1;

    for(uint n = 0; n < nMax; n++){
        m_R.push_back(zeros(tMax, tMax, tMax));
    }

    Boys boys(nMax);
    boys.evaluateBoysFunctions(alpha*dot(PQ,PQ));
    rowvec Fn = boys.getBoysFunctions();


    for(uint n = 0; n < nMax; n++){
        m_R.at(n)(0,0,0) = pow(-2*alpha,n)*Fn[n];
    }


    // p = previous
    // R(n,t,u,v) = (v-1) * R(n+1,t,u,v-2) + PQz * R(n+1,t,u,v-1)
    for(uint v = 1; v < vMax; v++){
        for(uint n = 0; n < nMax-v; n++){
            int t = 0.0; int u = 0.0;
            int vpp = v - 2;
            uint vp = v - 1;

            double R_t_u_vpp = 0.0;
            if(!(vpp < 0)){
                R_t_u_vpp = m_R.at(n+1)(t,u,vpp);
            }

            double R_t_u_vp = m_R.at(n+1)(t,u,vp);

            m_R.at(n)(t,u,v) = (v-1) * R_t_u_vpp + PQ(2) * R_t_u_vp;
        }
    }


    // p = previous
    // R(n,t,u,v) = (u-1) * R(n+1,t,u-2,v) + PQy * R(n+1,t,u-1,v)
    for(uint u = 1; u < uMax; u++){
        for(uint n = 0; n < nMax-u; n++){
            int t = 0.0; int v = 0.0;
            int upp = u - 2.0;
            uint up = u - 1.0;


            double R_t_upp_v = 0.0;
            if(!(upp < 0)){
                R_t_upp_v = m_R.at(n+1)(t,upp,v);
            }

            double R_t_up_v = m_R.at(n+1)(t,up,v);

            m_R.at(n)(t,u,v) = (u-1) * R_t_upp_v + PQ(1) * R_t_up_v;
        }
    }



    // p = previous
    // R(n,t,u,v) = (t-1) * R(n+1,t-2,u,v) + PQx * R(n+1,t-1,u,v)
    for(uint t = 1; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
                for(uint n = 0; n < nMax-t; n++){
                    int tpp = t - 2.0;
                    uint tp = t - 1.0;

                    double R_tpp_u_v = 0.0;
                    if(!(tpp < 0)){
                        R_tpp_u_v = m_R.at(n+1)(tpp,u,v);
                    }

                    double R_tp_u_v = m_R.at(n+1)(tp,u,v);

                    m_R.at(n)(t,u,v) = (t-1) * R_tpp_u_v + PQ(0) * R_tp_u_v;
                }
            }
        }
    }




}



void Integrator::setupE()
{
    int iAmax = m_maxAngularMomentum + 3;
    int iBmax = m_maxAngularMomentum + 3;
    int tmax  = 2*(m_maxAngularMomentum + 2) + 1;

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
            for(int t = 0; t < tmax; t++){

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

                m_E[cor](iA,iB,t) = 1.0 / (2*p) * E_iA_iBp_tp + PB(cor) * E_iA_iBp_t +  (t + 1)*E_iA_iBp_tn;
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

                    m_E[cor](iA,iB,t) = 1.0 / (2*p) * E_iAp_iB_tp + PA(cor) * E_iAp_iB_t +  (t + 1)*E_iAp_iB_tn;
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


    //        cout << m_E[0] << endl;
    //        cout << "-----------------------------------" <<endl;
    //        cout << m_E[1] << endl;
    //        cout << "-----------------------------------" <<endl;
    //        cout << m_E[2] << endl;

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


double Integrator::kineticIntegral(int cor, int iA, int iB) {
    double b = m_primitives[1]->exponent();

    double S_iA_iBnn = overlapIntegral(cor, iA, iB + 2);
    double S_iA_iB = overlapIntegral(cor, iA, iB);
    double S_iA_iBpp;
    if(iB - 2 >= 0) {
        S_iA_iBpp= overlapIntegral(cor, iA, iB - 2);
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


double Integrator::nuclearAttractionIntegral(int iA, int jA, int kA, int iB, int jB, int kB)
{


    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;

    double a = m_primitives[0]->exponent();
    double b = m_primitives[1]->exponent();

    double p = a + b;
    rowvec P  = (a*A + b*B)/p;
    rowvec PC = P - C;

    setupR(PC,p);

    double result = 0.0;

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
                result += m_E[0](iA, iB, t) * m_E[1](jA, jB, u) * m_E[2](kA, kB, v) * m_R.at(0)(t,u,v);
            }
        }
    }

    return 2 * result * M_PI / p;
}



double Integrator::electronRepulsionIntegral(int iA, int jA, int kA, int iB, int jB, int kB,
                                             int iC, int jC, int kC, int iD, int jD, int kD)
{


    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;
    const rowvec &C = m_corePositionC;
    const rowvec &D = m_corePositionD;

    double a = m_primitives[0]->exponent();
    double b = m_primitives[1]->exponent();
    double c = m_primitives[2]->exponent();
    double d = m_primitives[3]->exponent();

    double p = a + b;
    double q = c + d;
    rowvec P  = (a*A + b*B)/p;
    rowvec Q  = (c*C + d*D)/q;

    double alpha = p*q/(p+q);
    rowvec PQ = P - Q;

    setupR(PQ,alpha);

    double result = 0.0;
    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;
    uint kMax = iC + iD + 1;
    uint lMax = jC + jD + 1;
    uint mMax = kC + kD + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){

                double Etuv= m_E[0](iA, iB, t) * m_E[1](jA, jB, u) * m_E[2](kA, kB, v);

                for(uint k = 0; k < kMax; k++){
                    for(uint l = 0; l < lMax; l++){
                        for(uint m = 0; m < mMax; m++){

                            double Eklm= m_E[0](iC, iD, k) * m_E[1](jC, jD, l) * m_E[2](kC, kD, m);
                            result += Etuv * Eklm * m_R.at(0)(t+k,u+l,v+m) * pow(-1, k+l+m);
                        }
                    }
                }
            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q));

    return result;

}




















//    for(uint tuvSum = 0; tuvSum < nMax; tuvSum++){
//        for(uint n = 0; n < nMax-1; n++){

//            for(uint t = 0; t < tMax-1; t++){
//                for(uint u = 0; u < tMax-1; u++){
//                    for(uint v = 0; v < tMax-1; v++){

//                        if(t+u+v == tuvSum){
//                            int tp = t - 1;
//                            int up = u - 1;
//                            int vp = v - 1;

//                            double R_tp_u_v = 0.0;
//                            if(!(tp < 0)){
//                                R_tp_u_v = m_R.at(n+1)(tp,u,v);
//                            }

//                            double R_t_up_v = 0;
//                            if(!(up < 0)) {
//                                R_t_up_v = m_R.at(n+1)(t,up,v);
//                            }

//                            double R_t_u_vp = 0;
//                            if(!(vp < 0)) {
//                                R_t_u_vp = m_R.at(n+1)(t,u,vp);
//                            }

//                            double R_t_u_v = m_R.at(n+1)(t,u,v);

//                            m_R.at(n)(t+1,u,v) = t * R_tp_u_v + PQ(0) * R_t_u_v;
//                            m_R.at(n)(t,u+1,v) = u * R_t_up_v + PQ(1) * R_t_u_v;
//                            m_R.at(n)(t,u,v+1) = v * R_t_u_vp + PQ(2) * R_t_u_v;

////                            cout << n << t+1 << u << v << endl;
////                            cout << n << t << u+1 << v << endl;
////                            cout << n << t << u << v+1 << endl;

////                            cout <<  m_R.at(n)(t+1,u,v) << m_R.at(n)(t,u+1,v) << m_R.at(n)(t,u,v+1)<<endl;

//                        }
//                    }
//                }
//            }
////            cout << "----" <<endl;
////            sleep(3);


//        }
//    }


//    cout  << m_R.at(0) <<endl;





