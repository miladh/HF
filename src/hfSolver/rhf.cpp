#include "rhf.h"

using namespace hf;

RHF::RHF(ElectronicSystem *system):
    HFsolver(system),
    m_F(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_C(ones(m_nBasisFunctions,m_nBasisFunctions)),
    m_P(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_fockEnergy(zeros(m_nBasisFunctions))

{
}

void RHF::advance()
{
    vec fockEnergyOld;
    double stdDeviation= 1.0;
    int maxNumOfIteration = 200;
    m_iteration = 0;

    while (stdDeviation  > HFSOLVERTOLERANCE){
        if(m_rank == 0){
            cout << "Convergence rate: "
                 << 100.0*HFSOLVERTOLERANCE/stdDeviation
                 << setprecision(2)
                 <<" %" << endl;
        }
        fockEnergyOld = m_fockEnergy;
        solveSingle();
        stdDeviation = computeStdDeviation(m_fockEnergy, fockEnergyOld);
        updateFockMatrix();

        m_iteration+=1;
        if(m_iteration > maxNumOfIteration){
            if(m_rank == 0){
                cerr << "Energy has not converged! " << endl;
            }
            m_energy = 0.0;
            break;
        }
    }

}
void RHF::solveSingle()
{
    vec eigVal; mat eigVec, V;
    eig_sym(eigVal, eigVec, m_S);
    V = eigVec*diagmat(1.0/sqrt(eigVal));


//    DIISprocedure();

    eig_sym(eigVal, eigVec, V.t() * m_F * V);

    m_C = V*eigVec;

    normalize(m_C, m_nElectrons/2);

    m_P =  2.0 * m_C.cols(0, m_nElectrons/2.0-1) * m_C.cols(0, m_nElectrons/2.0-1).t();
    m_fockEnergy = eigVal;
}

void RHF::DIISprocedure()
{
    m_errors.push_back(m_F*m_P*m_S - m_S*m_P*m_F);
    m_fockMatrices.push_back(m_F);


    if(m_errors.size() > 20){
        m_errors.erase(m_errors.begin());
        m_fockMatrices.erase(m_fockMatrices.begin());
        mat A = zeros(m_errors.size()+1, m_errors.size()+1);

        for(uint i = 0; i < m_errors.size(); i++){
            A(i, m_errors.size()) = -1;
            A(m_errors.size(), i) = -1;

            for(uint j = 0; j < m_errors.size(); j++){
                A(i,j) = trace(m_errors.at(i) * m_errors.at(j));
            }
        }

        vec b = zeros(A.n_rows);
        b(b.n_elem -1) = -1.0;
        b = inv(A) * b;

        for(uint i = 0; i < b.n_elem - 1; i++){
            m_F += b(i) * m_fockMatrices.at(i);
        }
    }
}


void RHF::updateFockMatrix()
{
    for (int p = 0; p < m_nBasisFunctions; p++){
        for (int q = 0; q < m_nBasisFunctions; q++){

            m_F(p,q) = m_h(p,q);

            for (int r = 0; r < m_nBasisFunctions; r++){
                for (int s = 0; s < m_nBasisFunctions; s++){
                    m_F(p,q) += 0.5 * m_P(r,s) * (2.0 * m_Q(p,r)(q,s) - m_Q(p,r)(s,q));
                }
            }
        }
    }
}

void RHF::calculateEnergy()
{
    m_energy = 0.5 * accu(m_P % (m_h + m_F)) + m_system->nuclearPotential();
}

field<const mat *> RHF::fockMatrix()
{
    updateFockMatrix();
    field<const mat *> fockMatrices(1,1);
    fockMatrices(0,0) = &m_F;
    return fockMatrices;
}

field<const mat *> RHF::densityMatrix() const
{
    field<const mat *> densityMatrices(1,1);
    densityMatrices(0,0) = &m_P;
    return densityMatrices;
}

field<const mat *> RHF::expansionCoefficients() const
{
    field<const mat *> coefficientsMatrices(1,1);
    coefficientsMatrices(0,0) = &m_C;
    return coefficientsMatrices;
}
