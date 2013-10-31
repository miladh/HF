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
    for(int step=0; step < m_nSteps; step++){

        setupTwoParticleMatrix();
        m_F = m_h + m_G;

        vec s; mat U;
        eig_sym(s, U, m_S);

        mat V = U*diagmat(1.0/sqrt(s));
        m_F = V.t() * m_F * V;


        cout << m_S <<endl;
        vec eps;
        mat Cmat;
        eig_sym(eps, Cmat, m_F);

        m_C = V*Cmat.col(0);
        normalize();

        double Eg=0.0;

        for(uint p=0; p < m_C.n_elem; p++){
            for(uint q=0; q < m_C.n_elem; q++){
                Eg += m_C(p)*m_C(q)*m_h(p,q);
            }
        }

        Eg = 2*Eg;
        for(uint p=0; p < m_C.n_elem; p++){
            for(uint r=0; r < m_C.n_elem; r++){
                for(uint q=0; q< m_C.n_elem; q++){
                    for(uint s=0; s < m_C.n_elem; s++){
                        Eg +=m_Q[p][r][q][s]*m_C(p)*m_C(q)*m_C(r)*m_C(s);
                    }
                }
            }
        }

        cout <<"Energy: " << Eg <<" step: " << step << endl;

    }



}
void HFsolver::setupTwoParticleMatrix()
{

    for(uint a=0; a < m_C.n_elem; a++){
        for(uint b=0; b < m_C.n_elem; b++){
            for(uint c=0; c < m_C.n_elem; c++){
                for(uint d=0; d < m_C.n_elem; d++){
                    m_G(a,b) += m_Q[a][c][b][d]*m_C(c)*m_C(d);
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
