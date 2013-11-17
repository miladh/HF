#include "hfsolver.h"

HFsolver::HFsolver(System *system):
    m_system(system),
    m_nElectrons(system->getNumOfElectrons()),
    m_nOrbitals(system->getTotalNumOfBasisFunc()),
    m_S(zeros(m_nOrbitals,m_nOrbitals)),
    m_h(zeros(m_nOrbitals,m_nOrbitals)),
    m_F(zeros(m_nOrbitals,m_nOrbitals)),
    m_P(zeros(m_nOrbitals,m_nOrbitals)),
    m_C(ones(m_nOrbitals,m_nElectrons/2.0))

{
    m_Q.set_size(m_nOrbitals, m_nOrbitals);
    for(int i = 0; i < m_nOrbitals; i++){
        for(int j = 0; j < m_nOrbitals; j++){
            m_Q(i,j) = zeros(m_nOrbitals,m_nOrbitals);
        }
    }

}

void HFsolver::runSolver()
{
    double fockEnergyOld;
    double energyDiff = 1.0;
    int step = 0;
    int maxStep = 100;

    setupOneParticleMatrix();
    setupTwoParticleMatrix();

    while (energyDiff > HFSOLVERTOLERANCE){
        fockEnergyOld = m_fockEnergy;
        setupFockMatrix();
        solveSingle();
        energyDiff = fabs(fockEnergyOld - m_fockEnergy);
        step+=1;
        if(step > maxStep){
            cerr << "Energy has not converged! " << endl;
        }
    }

    calculateEnergy();
}

void HFsolver::solveSingle()
{
    vec eigVal;
    mat eigVec;
    eig_sym(eigVal, eigVec, m_S);

    mat V = eigVec*diagmat(1.0/sqrt(eigVal));

    m_F = V.t()*m_F*V;


    eig_sym(eigVal, eigVec, m_F);
    m_C = V*eigVec.cols(0, m_nElectrons/2.0-1);

    normalize();

    m_P = 2*m_C*m_C.t();

    m_fockEnergy = eigVal(0);
}

void HFsolver::calculateEnergy()
{
    m_energy = 0;

    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){
            m_energy += m_P(p, q)*m_h(p, q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    m_energy += 0.5*m_P(p,q)*m_P(s,r)*(m_Q(p,r)(q,s) - 0.5*m_Q(p,r)(s,q));
                }
            }
        }
    }
    m_energy += m_system->getNucleiPotential();

//    cout << "Energy: " << setprecision(10) << m_energy << endl;

}

void HFsolver::setupOneParticleMatrix()
{
    rowvec oneElectronIntegrals;

    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            oneElectronIntegrals = m_system->getOneParticleIntegral(p,q);
            m_S(p,q) = oneElectronIntegrals(0);
            m_h(p,q) = oneElectronIntegrals(1);
        }
    }

//    m_S = symmatu(m_S);
//    m_h = symmatu(m_h);

}


void HFsolver::setupTwoParticleMatrix()
{
    for(int p = 0; p < m_nOrbitals; p++){
        for(int r = 0; r < m_nOrbitals; r++){
            for(int q = 0; q < m_nOrbitals; q++){
                for(int s = 0; s < m_nOrbitals; s++){

                    m_Q(p,r)(q,s) = m_system->getTwoParticleIntegral(p,q,r,s);
                }
            }
        }
    }

}

void HFsolver::setupFockMatrix()
{
    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){

            m_F(p,q) = m_h(p,q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    m_F(p,q) += 0.5 * m_P(s,r) * (2 * m_Q(p,r)(q,s) - m_Q(p,r)(s,q));
                }
            }
        }
    }
}

void HFsolver::normalize()
{
    double norm;
    for (int i = 0; i < m_nElectrons/2; i++){
        norm = dot(m_C.col(i), m_S*m_C.col(i));
        m_C.col(i) = m_C.col(i)/sqrt(norm);
    }
}

double HFsolver::getEnergy() const
{
    return m_energy;
}

field<mat> HFsolver::getQmatrix(){

    return m_Q;
}
mat HFsolver::gethmatrix(){

    return m_h;
}

mat HFsolver::getSmatrix(){

    return m_S;
}
mat HFsolver::getC() const
{
    return m_C;
}








