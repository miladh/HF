#include "hfsolver.h"

HFsolver::HFsolver(System system, int nOrbitals, int nSteps):
    m_system(system),
    m_nSteps(nSteps),
    m_F(zeros(nOrbitals,nOrbitals)),
    m_S(zeros(nOrbitals,nOrbitals)),
    m_G(zeros(nOrbitals,nOrbitals)),
    m_C(ones(nOrbitals,nOrbitals)),
    m_P(zeros(nOrbitals,nOrbitals))

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

        setupDensityMatrix();
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

        for(uint a=0; a < m_P.n_rows; a++){
            for(uint b=0; b < m_P.n_cols; b++){
                Eg += m_P(a,b)*m_h(a,b);
            }
        }

        for(uint a=0; a < m_P.n_rows; a++){
            for(uint b=0; b < m_P.n_cols; b++){
                for(uint c=0; c< m_P.n_rows; c++){
                    for(uint d=0; d < m_P.n_cols; d++){
                        Eg +=m_Q(a,c)(b,d)*m_P(a,b)*m_P(c,d);
                    }
                }
            }
        }

        Eg *=0.5;

        cout.precision(8);
        cout <<"Energy: " << Eg <<" step: " << step << endl;

    }



}
void HFsolver::setupTwoParticleMatrix()
{
    for(uint a=0; a < m_G.n_rows; a++){
        for(uint b=0; b < m_G.n_cols; b++){

            for(uint c=0; c < m_P.n_rows; c++){
                for(uint d=0; d < m_P.n_cols; d++){
                    m_G(a,b) += m_Q(a,b)(c,d)*m_P(c,d);

                }
            }
        }
    }

}

void HFsolver::setupDensityMatrix()
{
    for(int p = 0; p < m_P.n_rows; p++){
        for(int q = 0; q < m_P.n_cols; q++){

            for(int k = 0; k < 2*0.5; k++){

                m_P(p,q) += 2*m_C(p,k)*m_C(q,k);

            }

        }

    }
}

void HFsolver::normalize(){
    for(int k=0; k < m_C.n_cols; k++){
        double normFactor= 0.0;

        for(uint i= 0; i < m_C.n_rows; i++){
            for(uint j= 0; j < m_C.n_rows; j++){
                normFactor += m_C(i,k)*m_S(i,j)*m_C(j,k);
            }
        }

        m_C.col(k) /= sqrt(normFactor);
    }
}
