#include "hfsolver.h"

HFsolver::HFsolver(System system):
    m_system(system),
    m_nElectrons(system.getNumOfElectrons()),
    m_nOrbitals(system.getTotalNumOfBasisFunc()),
    m_F(zeros(m_nOrbitals,m_nOrbitals)),
    m_S(zeros(m_nOrbitals,m_nOrbitals)),
    m_P(zeros(m_nOrbitals,m_nOrbitals)),
    m_C(zeros(m_nOrbitals,m_nElectrons/2.0))

{
    m_S = m_system.getOverlapMatrix();
    m_h = m_system.getOneParticleMatrix();
    m_Q = m_system.getTwoParticleMatrix();
    normalize();

    m_fockEnergy = 1.0E6;
    m_energy = 1.0E6;
    m_toler = 1.0E-6;
}

void HFsolver::runSolver()
{

    double m_fockEnergyOld;
    double m_energyDiff = 1.0;

    m_C = zeros<mat>(m_nOrbitals, m_nElectrons/2.0);

    // Iterate until the fock m_energy has converged
    while (m_energyDiff > m_toler){
        m_fockEnergyOld = m_fockEnergy;
        setupFockMatrix();
        solveSingle();
        m_energyDiff = fabs(m_fockEnergyOld - m_fockEnergy);

        // Calculate m_energy (not equal to Fock m_energy)
        m_energy = 0;

        for (int a = 0; a < m_nOrbitals; a++){
            for (int b = 0; b < m_nOrbitals; b++){
                m_energy += m_P(a, b)*m_h(a, b);

                for (int c = 0; c < m_nOrbitals; c++){
                    for (int d = 0; d < m_nOrbitals; d++){
                        m_energy += 0.5*m_P(a,b)*m_P(c,d)*(m_Q(a,b)(c,d) - 0.5*m_Q(a,b)(d,c));
                    }
                }
            }
        }

        cout << "Energy: " << setprecision(10) << m_energy << endl;

    }


}


void HFsolver::setupFockMatrix()
{

    for (int a = 0; a < m_nOrbitals; a++){
        for (int b = 0; b < m_nOrbitals; b++){

            // One-electron integrals
            m_F(a,b) = m_h(a,b);

            // Add two-electron integrals
            for (int c = 0; c < m_nOrbitals; c++){
                for (int d = 0; d < m_nOrbitals; d++){
                    m_F(a,b) += 0.5*m_P(c,d)*(2*m_Q(a,b)(c,d) - m_Q(a,b)(d,c));
                }
            }
        }
    }
}

void HFsolver::solveSingle()
{
    vec eigVal;
    mat eigVec;
    eig_sym(eigVal, eigVec, m_S);

    mat V = eigVec*diagmat(1.0/sqrt(eigVal));

    m_F = V.t()*m_F*V;


    eig_sym(eigVal, eigVec, m_F);
    m_C = V*eigVec.cols(0, m_nElectrons/2-1);

    // Normalize vector C
    normalize();

    // Compute density matrix
    m_P = 2*m_C*m_C.t();

    m_fockEnergy = eigVal(0);
}


void HFsolver::normalize()
{
    double norm;
    for (int i = 0; i < m_nElectrons/2; i++){
        norm = dot(m_C.col(i), m_S*m_C.col(i));
        m_C.col(i) = m_C.col(i)/sqrt(norm);
    }

}



