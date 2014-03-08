#include "uhf.h"

using namespace hf;

UHF::UHF(System *system, const int &rank, const int &nProcs):
    HFsolver(system, rank, nProcs),
    m_Fu(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Cu(ones(m_nBasisFunctions,m_nBasisFunctions)),
    m_Pu(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Fd(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Cd(ones(m_nBasisFunctions,m_nBasisFunctions)),
    m_Pd(zeros(m_nBasisFunctions,m_nBasisFunctions))
{
    m_nSpinUpElectrons   = 5;
    m_nSpinDownElectrons = 5;
    m_Pu(0,0) = 1.0;
}



void UHF::advance()
{
    double fockEnergyUOld, fockEnergyDOld;
    double energyDiff = 1.0;
    int step = 0;
    int maxStep = 100;

    while (energyDiff > HFSOLVERTOLERANCE){
        fockEnergyUOld = m_fockEnergyU;
        fockEnergyDOld = m_fockEnergyD;
        solveSingle();
        energyDiff = max(fabs(fockEnergyUOld - m_fockEnergyU),fabs(fockEnergyDOld - m_fockEnergyD));
        updateFockMatrix();

        step+=1;
        if(step > maxStep){
            cerr << "Energy has not converged! " << endl;
            exit(1);
        }
    }

}


void UHF::solveSingle()
{
    vec eigVal;
    mat eigVec, V;

    eig_sym(eigVal, eigVec, m_S);
    V = eigVec*diagmat(1.0/sqrt(eigVal));

    eig_sym(eigVal, eigVec, V.t()*m_Fu*V);
    m_Cu = V*eigVec;
    m_fockEnergyU = eigVal(0);

    eig_sym(eigVal, eigVec, V.t()*m_Fd*V);
    m_Cd = V*eigVec;
    m_fockEnergyD = eigVal(0);


    normalize();

    m_Pu = 0.5 * m_Pu + 0.5*m_Cu.cols(0, m_nSpinUpElectrons-1) * m_Cu.cols(0, m_nSpinUpElectrons-1).t();
    m_Pd = 0.5 * m_Pd + 0.5*m_Cd.cols(0, m_nSpinDownElectrons-1) * m_Cd.cols(0, m_nSpinDownElectrons-1).t();
}

void UHF::updateFockMatrix()
{
    for (int p = 0; p < m_nBasisFunctions; p++){
        for (int q = 0; q < m_nBasisFunctions; q++){

            m_Fu(p,q) = m_h(p,q);
            m_Fd(p,q) = m_h(p,q);

            for (int r = 0; r < m_nBasisFunctions; r++){
                for (int s = 0; s < m_nBasisFunctions; s++){

                    m_Fu(p,q) += m_Pu(s,r) * (m_Q(p,r)(q,s) - m_Q(p,r)(s,q)) + m_Pd(s,r) * m_Q(p,r)(q,s);
                    m_Fd(p,q) += m_Pd(s,r) * (m_Q(p,r)(q,s) - m_Q(p,r)(s,q)) + m_Pu(s,r) * m_Q(p,r)(q,s);
                }
            }
        }
    }
}

void UHF::calculateEnergy()
{
    m_energy = 0.5 * accu( (m_Pu + m_Pd) % m_h + m_Fu % m_Pu + m_Fd % m_Pd)
               + m_system->getNucleiPotential();
}


void UHF::normalize()
{
    double norm;
    for (int i = 0; i < m_nSpinUpElectrons; i++){
        norm = dot(m_Cu.col(i), m_S * m_Cu.col(i));
        m_Cu.col(i) = m_Cu.col(i)/sqrt(norm);
    }

    for (int i = 0; i < m_nSpinDownElectrons; i++){
        norm = dot(m_Cd.col(i), m_S * m_Cd.col(i));
        m_Cd.col(i) = m_Cd.col(i)/sqrt(norm);
    }
}


void UHF::calculateDensity()
{

    cout << "---Calculating density---" << endl;

    vec x = linspace(-10, 10, m_nProcs * 10);
    vec y = linspace(-10, 10, m_nProcs * 10);
    vec z = linspace(-10, 10, m_nProcs * 10);
    double dr = (x(1) - x(0)) * (y(1) - y(0)) * (z(1) - z(0));

    m_density = zeros(x.n_elem, y.n_elem, z.n_elem);
    double sumDensity =0;


    int xElements = x.n_elem/m_nProcs;
    int yElements = y.n_elem/m_nProcs;
    int zElements = z.n_elem/m_nProcs;

    int xMin = m_rank % m_nProcs * xElements;
    int yMin = m_rank / m_nProcs * yElements;
    int zMin = m_rank / m_nProcs * zElements;

    int xMax = xMin + xElements;
    int yMax = m_nProcs * yElements;
    int zMax = m_nProcs *  zElements;

    for(int i = xMin; i < xMax; i++) {
        for(int j = yMin; j < yMax; j++) {
            for(int k = zMin; k < zMax; k++) {

                for(int p = 0; p < m_nBasisFunctions; p++){
                    double innerProduct = m_system->gaussianProduct(p, p, x(i), y(j), z(k));
                    sumDensity += (m_Pu(p,p) + m_Pu(p,p)) * innerProduct * dr;
                    m_density(j,i,k) += (m_Pu(p,p) + m_Pu(p,p)) * innerProduct ;
                    m_density(j,i,k) += (m_Pu(p,p) + m_Pu(p,p)) * innerProduct ;

                    for(int q = p+1; q < m_nBasisFunctions; q++){
                        innerProduct = m_system->gaussianProduct(p, q, x(i), y(j), z(k));
                        sumDensity += 2.0 * (m_Pu(p,p) + m_Pu(p,p)) * innerProduct * dr;
                        m_density(j,i,k) += 2.0 * (m_Pu(p,p) + m_Pu(p,p)) * innerProduct ;

                    }
                }

            }
        }
    }

    cout << "density sum: " << sumDensity << endl;
    densityOutput(x.min(),x.max(),y.min(),y.max(),z.min(),z.max());

}
