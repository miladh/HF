#include "geometricalderivative.h"

using namespace hf;

GeometricalDerivative::GeometricalDerivative(System *system, RHF *solver):
    m_system(system),
    m_solver(solver),
    m_nBasisFunctions(system->getTotalNumOfBasisFunc())
{
    //initialize:
    m_dh.set_size(m_nBasisFunctions,m_nBasisFunctions);
    m_dS.set_size(m_nBasisFunctions,m_nBasisFunctions);
    m_dQ.set_size(m_nBasisFunctions, m_nBasisFunctions);
    for(int i = 0; i < m_nBasisFunctions; i++ ){
        for(int j = 0; j < m_nBasisFunctions; j++ ){
            m_dh(i,j)  = zeros<rowvec>(3);
            m_dS(i,j)  = zeros<rowvec>(3);
            m_dQ(i,j).set_size(m_nBasisFunctions,m_nBasisFunctions);
        }
    }

    for(int i = 0; i < m_nBasisFunctions; i++ ){
        for(int j = 0; j < m_nBasisFunctions; j++ ){
            for(int k = 0; k < m_nBasisFunctions; k++ ){
                for(int l = 0; l < m_nBasisFunctions; l++ ){
                    m_dQ(i,j)(k,l) = zeros<rowvec>(3);
                }
            }
        }
    }


}


rowvec GeometricalDerivative::energyGradient(const int core)
{
    m_differentiationCore = core;
    setupDerivativeMatrices();
    return calculateEnergyGradient();

}



void GeometricalDerivative::setupDerivativeMatrices()
{

    mat diffOneParticleIntegral;
    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int q = 0; q < m_nBasisFunctions; q++){
            diffOneParticleIntegral = m_system->getOneParticleDerivative(p,q,m_differentiationCore);
            m_dS(p,q) = diffOneParticleIntegral.row(0);
            m_dh(p,q) = diffOneParticleIntegral.row(1);
        }
    }


    //Set up the dQ array:
    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int r = 0; r < m_nBasisFunctions; r++){
            for(int q = 0; q < m_nBasisFunctions; q++){
                for(int s = 0; s < m_nBasisFunctions; s++){
                    m_dQ(p,r)(q,s) = m_system->getTwoParticleIntegralDerivative(p,q,r,s,m_differentiationCore);

                }
            }
        }
    }

}

rowvec GeometricalDerivative::calculateEnergyGradient()
{
    rowvec dE  = {0,0,0};
    mat F = m_solver->getFockMatrix();
    mat P = m_solver->getDensityMatrix();

    for (int p = 0; p < m_nBasisFunctions; p++){
        for (int q = 0; q < m_nBasisFunctions; q++){
            dE += P(p, q)*m_dh(p, q);

            for (int r = 0; r < m_nBasisFunctions; r++){
                for (int s = 0; s < m_nBasisFunctions; s++){
                    dE += 0.5*P(p,q)*P(s,r)*(m_dQ(p,r)(q,s) - 0.5*m_dQ(p,r)(s,q));
                }
            }
        }
    }


    mat dSx, dSy,dSz;
    dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
    dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
    dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);
    for(int i = 0; i < m_nBasisFunctions; i++){
        for(int j = 0; j < m_nBasisFunctions; j++){
            dSx(i,j) = m_dS(i, j)(0);
            dSy(i,j) = m_dS(i, j)(1);
            dSz(i,j) = m_dS(i, j)(2);
        }
    }

    dE(0) -= 0.5 * trace(P*dSx*P*F);
    dE(1) -= 0.5 * trace(P*dSy*P*F);
    dE(2) -= 0.5 * trace(P*dSz*P*F);

    //    Nuclear repulsion term
    dE  +=m_system->getNucleiPotential_derivative(m_differentiationCore);

    return dE;
}

