#include "integrator.h"


Integrator::Integrator():
    m_corePositionA(zeros<rowvec>(3)),
    m_corePositionB(zeros<rowvec>(3))
{
    //    setMaxAngularMomentum(0);
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


uint Integrator::maxAngularMomentum() const
{
    return m_maxAngularMomentum;
}

void Integrator::setMaxAngularMomentum(const uint &maxAngularMomentum)
{
    m_maxAngularMomentum = maxAngularMomentum;

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

    m_boys = new Boys(nMax_ee - 1);
    m_E.set_size(3);
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


void Integrator::setupR(const rowvec &PQ, const double &alpha, field<cube> &R){


    int tMax  = R(0).n_rows   - 1;
    int uMax  = R(0).n_cols   - 1;
    int vMax  = R(0).n_slices - 1;
    int nMax  = R.n_elem      - 1;

    m_boys->evaluateBoysFunctions(alpha*dot(PQ,PQ));

    for(int n = 0; n < nMax+1; n++){
        R(n)(0,0,0) = pow(-2*alpha,n)*m_boys->getBoysFunctions(n);
    }

    // p = previous
    // R(n,t,u,v+1) = v * R(n+1,t,u,v-1) + PQz * R(n+1,t,u,v)
    for(int v = 0; v < vMax; v++){ //not including vMax, since we have v+1 in formula
        for(int n = 0; n < nMax-v; n++){//not including nMax, since we have n+1 in formula
            int t = 0.0; int u = 0.0;
            int vp = v - 1.0;

            double R_t_u_vp = 0.0;
            if(!(vp < 0)){
                R_t_u_vp = R(n+1)(t,u,vp);
            }

            double R_t_u_v = R(n+1)(t,u,v);

            R(n)(t,u,v+1) = v * R_t_u_vp + PQ(2) * R_t_u_v;
        }
    }


    // p = previous
    // R(n,t,u+1,v) = u * R(n+1,t,u-1,v) + PQy * R(n+1,t,u,v)
    for(int v = 0; v < vMax+1; v++){ //including vMax, since we have v in formula (not v+1)
        for(int u = 0; u < uMax - v; u++){ //not including uMax, since we have u+1 in formula
            for(int n = 0; n < nMax - u - v ; n++){//not including nMax, since we have n+1 in formula
                int t = 0.0;
                int up = u - 1.0;

                double R_t_up_v = 0.0;
                if(!(up < 0)){
                    R_t_up_v = R(n+1)(t,up,v);
                }

                double R_t_u_v = R(n+1)(t,u,v);
                R(n)(t,u+1,v) = u * R_t_up_v + PQ(1) * R_t_u_v;
            }
        }
    }


    // p = previous
    // R(n,t+1,u,v) = t * R(n+1,t-1,u,v) + PQx * R(n+1,t,u,v)
    for(int u = 0; u < uMax+1; u++){//including uMax, since we have u in formula (not u+1)
        for(int v = 0; v < vMax+1; v++){//including vMax, since we have v in formula (not v+1)
            for(int t = 0; t < tMax - u - v; t++){ //not including tMax, since we have t+1 in formula
                for(int n = 0; n < nMax - t - u - v; n++){//not including nMax, since we have n+1 in formula
                    int tp = t - 1.0;

                    double R_tp_u_v = 0.0;
                    if(!(tp < 0)){
                        R_tp_u_v = R(n+1)(tp,u,v);
                    }

                    double R_t_u_v = R(n+1)(t,u,v);

                    R(n)(t+1,u,v) = t * R_tp_u_v + PQ(0) * R_t_u_v;
                }
            }
        }
    }
}


void Integrator::setupE()
{
    uint iAmax = m_maxAngularMomentum + 3;
    uint iBmax = m_maxAngularMomentum + 3;
    uint tmax  = 2*(m_maxAngularMomentum + 2) + 1;

    for(uint cor = 0; cor < 3; cor++){
        m_E[cor] = zeros(iAmax, iBmax, tmax);
    }

    const rowvec &A = m_corePositionA;
    const rowvec &B = m_corePositionB;

    double a = m_exponentA;
    double b = m_exponentB;

    double p = a + b;
    double mu = a * b / p;

    rowvec P  = (a*A + b*B)/p;
    rowvec AB = A - B;
    rowvec PA = P - A;
    rowvec PB = P - B;
    rowvec Kab = exp(-mu*AB%AB);


    for(uint cor = 0; cor < 3; cor++){
        m_E[cor](0,0,0) = Kab(cor);
    }

    for(uint cor=0; cor < 3; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
        for(uint iB = 1; iB < iBmax; iB++){
            for(uint t = 0; t < tmax; t++){

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
        for(uint iA = 1; iA < iAmax; iA++) {
            for(uint iB = 0; iB < iBmax; iB++) {
                for(uint t = 0; t < tmax; t++) {

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
    double a = m_exponentA;
    double b = m_exponentB;
    double p = a + b;
    return m_E[cor](iA,iB,0) * sqrt(M_PI / p);
}


double Integrator::overlapIntegral(int iA, int jA, int kA, int iB, int jB, int kB) {

    return overlapIntegral(0, iA, iB) * overlapIntegral(1, jA, jB) * overlapIntegral(2, kA, kB);
}


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

    setupR(PC,p, m_Ren);

    double result = 0.0;

    uint tMax = iA + iB + 1;
    uint uMax = jA + jB + 1;
    uint vMax = kA + kB + 1;

    for(uint t = 0; t < tMax; t++){
        for(uint u = 0; u < uMax; u++){
            for(uint v = 0; v < vMax; v++){
                result += m_E[0](iA, iB, t) * m_E[1](jA, jB, u) * m_E[2](kA, kB, v) * m_Ren(0)(t,u,v);
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

    setupR(PQ,alpha, m_Ree);


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
                setupE(A,B,a,b);


                double Etuv= m_E[0](iA, iB, t) * m_E[1](jA, jB, u) * m_E[2](kA, kB, v);
                setupE(C,D,c,d);


                for(int k = 0; k < kMax; k++){
                    for(int l = 0; l < lMax; l++){
                        for(int m = 0; m < mMax; m++){

                            double Eklm= m_E[0](iC, iD, k) * m_E[1](jC, jD, l) * m_E[2](kC, kD, m);
                            result += Eklm * m_Ree(0)(t+k,u+l,v+m) * (1 - 2* ((k+l+m)%2));
                        }
                    }
                }


                result *=Etuv;
            }
        }
    }

    result *= 2*pow(M_PI,2.5)/ (p*q*sqrt(p+q));

//    cout << A(0) << "  " << B(0) << " --- " <<C(0) << "  " << D(0) <<endl;
//    cout << a << "  " << b << " --- " << c << "  " << d <<endl;
//    cout << "  res: " << result <<endl;
//    cout << endl;
//    if (D(0)==0.5){
//        cerr << "WARNING!!!" << endl <<endl;

//            sleep(5);
//    }

    return result;

}





void Integrator::setupE(const rowvec &A , const rowvec &B, const double &a, const double&b)
{
    uint iAmax = m_maxAngularMomentum + 3;
    uint iBmax = m_maxAngularMomentum + 3;
    uint tmax  = 2*(m_maxAngularMomentum + 2) + 1;

    for(uint cor = 0; cor < 3; cor++){
        m_E[cor] = zeros(iAmax, iBmax, tmax);
    }



    double p = a + b;
    double mu = a * b / p;

    rowvec P  = (a*A + b*B)/p;
    rowvec AB = A - B;
    rowvec PA = P - A;
    rowvec PB = P - B;
    rowvec Kab = exp(-mu*AB%AB);


    for(uint cor = 0; cor < 3; cor++){
        m_E[cor](0,0,0) = Kab(cor);
    }

    for(uint cor=0; cor < 3; cor++){ //Loop for x,y,z


        // p = previous, n = next
        // E(t,i,j) = 1 / (2*p) * E(t-1,i,j-1) + XPA * E(t,i,j-1) + (t + 1)*E(t+1,i,j-1)
        for(uint iB = 1; iB < iBmax; iB++){
            for(uint t = 0; t < tmax; t++){

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
        for(uint iA = 1; iA < iAmax; iA++) {
            for(uint iB = 0; iB < iBmax; iB++) {
                for(uint t = 0; t < tmax; t++) {

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














