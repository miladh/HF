#include "hfsolver.h"

HFsolver::HFsolver(System system, int nOrbitals, int nSteps):
    m_system(system),
    m_nSteps(nSteps),
    m_F(zeros(nOrbitals,nOrbitals)),
    m_S(zeros(nOrbitals,nOrbitals)),
    m_G(zeros(nOrbitals,nOrbitals)),
    m_C(ones(nOrbitals))
{
    m_S = m_system.getOverlapMatrix();
    m_h = m_system.getOneParticleMatrix();
    m_Q = m_system.getTwoParticleMatrix();
    normalize();
}

void HFsolver::runSolver()
{

//    for(int i = 0; i <m_C.n_elem; i++ ){
//        for(int j = 0; j <m_C.n_elem; j++ ){
//            for(int k = 0; k <m_C.n_elem; k++ ){
//                for(int l = 0; l <m_C.n_elem; l++ ){

//                    cout << m_Q(i,j)(k,l) <<endl;
//                }
//            }
//        }

//    }

//    sleep(7);


    for(int step=0; step < m_nSteps; step++){
        m_G = 0*m_G;
        setupTwoParticleMatrix();

        m_F = m_h + m_G;

        vec s; mat U;
        eig_sym(s, U, m_S);

        mat V = U*diagmat(1.0/sqrt(s));
        m_F = V.t() * m_F * V;


        vec eps;
        mat Cmat;
        eig_sym(eps, Cmat, m_F);

        m_C = V*Cmat.col(0);
        normalize();

        double Eg=0.0;

        for(uint a=0; a < m_C.n_elem; a++){
            for(uint b=0; b < m_C.n_elem; b++){
                Eg += m_C(a)*m_C(b)*m_h(a,b);
            }
        }

        Eg = 2*Eg;
        for(uint a=0; a < m_C.n_elem; a++){
            for(uint b=0; b < m_C.n_elem; b++){
                for(uint c=0; c< m_C.n_elem; c++){
                    for(uint d=0; d < m_C.n_elem; d++){
                        Eg +=m_Q(a,c)(b,d)*m_C(a)*m_C(b)*m_C(c)*m_C(d);
                    }
                }
            }
        }
         cout.precision(8);
        cout <<"Energy: " << Eg <<" step: " << step << endl;

    }



}
void HFsolver::setupTwoParticleMatrix()
{
    for(uint a=0; a < m_C.n_elem; a++){
        for(uint b=0; b < m_C.n_elem; b++){

            for(uint c=0; c < m_C.n_elem; c++){
                for(uint d=0; d < m_C.n_elem; d++){
                    m_G(a,b) += m_Q(a,b)(c,d)*m_C(c)*m_C(d);

                }
            }
        }
    }

}

void HFsolver::normalize(){
    double normFactor= 0.0;

    for(uint i= 0; i < m_C.n_elem; i++){
        for(uint j= 0; j < m_C.n_elem; j++){
            normFactor += m_C(i)*m_S(i,j)*m_C(j);
        }
    }


    m_C /= sqrt(normFactor);

}
