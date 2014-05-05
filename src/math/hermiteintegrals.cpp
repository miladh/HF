#include "hermiteintegrals.h"
using namespace hf;
HermiteIntegrals::HermiteIntegrals(const int highestOrder):
    m_boys(new Boys(highestOrder))
{
    m_R.set_size(highestOrder);
    for(int n = 0; n < highestOrder; n++){
        m_R(n) = zeros(highestOrder, highestOrder, highestOrder);
    }
}


double HermiteIntegrals::R(const int n,  const int t, const int u, const int v)const
{
    return m_R(n)(t,u,v);
}

void HermiteIntegrals::updateR(const rowvec &PQ, const double &alpha,
                              const int tMax, const int uMax, const int vMax){


    int nMax  = std::min(tMax + uMax + vMax, signed(m_R.n_elem - 1));
    m_boys->evaluateBoysFunctions(alpha*dot(PQ,PQ));

    for(int n = 0; n < nMax+1; n++){
        m_R(n)(0,0,0) = std::pow(-2*alpha,n)*m_boys->getBoysFunctions(n);
    }

    // p = previous
    // R(n,t,u,v+1) = v * R(n+1,t,u,v-1) + PQz * R(n+1,t,u,v)
    for(int v = 0; v < vMax; v++){ //not including vMax, since we have v+1 in formula
        for(int n = 0; n < nMax-v; n++){//not including nMax, since we have n+1 in formula
            int t = 0.0; int u = 0.0;
            int vp = v - 1.0;

            double R_t_u_vp = 0.0;
            if(!(vp < 0)){
                R_t_u_vp = m_R(n+1)(t,u,vp);
            }

            double R_t_u_v = m_R(n+1)(t,u,v);

            m_R(n)(t,u,v+1) = v * R_t_u_vp + PQ(2) * R_t_u_v;
        }
    }


    // p = previous
    // R(n,t,u+1,v) = u * R(n+1,t,u-1,v) + PQy * R(n+1,t,u,v)
    for(int v = 0; v < vMax+1; v++){ //including vMax, since we have v in formula (not v+1)
        for(int u = 0; u < uMax; u++){ //not including uMax, since we have u+1 in formula
            for(int n = 0; n < nMax - u - v ; n++){//not including nMax, since we have n+1 in formula
                int t = 0.0;
                int up = u - 1.0;

                double R_t_up_v = 0.0;
                if(!(up < 0)){
                    R_t_up_v = m_R(n+1)(t,up,v);
                }

                double R_t_u_v = m_R(n+1)(t,u,v);
                m_R(n)(t,u+1,v) = u * R_t_up_v + PQ(1) * R_t_u_v;
            }
        }
    }


    // p = previous
    // R(n,t+1,u,v) = t * R(n+1,t-1,u,v) + PQx * R(n+1,t,u,v)
    for(int u = 0; u < uMax+1; u++){//including uMax, since we have u in formula (not u+1)
        for(int v = 0; v < vMax+1; v++){//including vMax, since we have v in formula (not v+1)
            for(int t = 0; t < tMax ; t++){ //not including tMax, since we have t+1 in formula
                for(int n = 0; n < nMax - t - u - v; n++){//not including nMax, since we have n+1 in formula
                    int tp = t - 1.0;

                    double R_tp_u_v = 0.0;
                    if(!(tp < 0)){
                        R_tp_u_v = m_R(n+1)(tp,u,v);
                    }

                    double R_t_u_v = m_R(n+1)(t,u,v);

                    m_R(n)(t+1,u,v) = t * R_tp_u_v + PQ(0) * R_t_u_v;
                }
            }
        }
    }
}




