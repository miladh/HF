#include "uhf.h"

using namespace hf;

UHF::UHF(System *system, const int &rank, const int &nProcs):
    HFsolver(system, rank, nProcs),
    m_Fu(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Cu(ones(m_nBasisFunctions,m_nBasisFunctions)),
    m_Pu(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Fd(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_Cd(ones(m_nBasisFunctions,m_nBasisFunctions)),
    m_Pd(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_fockEnergyU(zeros(m_nBasisFunctions)),
    m_fockEnergyD(zeros(m_nBasisFunctions))
{
    m_Pu(0,0) = 0.1;
}


void UHF::advance()
{
    vec fockEnergyOld_U, fockEnergyOld_D;
    double stdDeviation_U = 1.0;
    double stdDeviation_D = 1.0;
    int maxNumOfIteration = 200;
    m_iteration = 0;

    while (stdDeviation_U > HFSOLVERTOLERANCE && stdDeviation_D > HFSOLVERTOLERANCE){
        fockEnergyOld_U = m_fockEnergyU;
        fockEnergyOld_D = m_fockEnergyD;
        solveSingle();
        updateFockMatrix();

        stdDeviation_U = computeStdDeviation(m_fockEnergyU, fockEnergyOld_U);
        stdDeviation_D = computeStdDeviation(m_fockEnergyD, fockEnergyOld_D);

        m_iteration+=1;
        if(m_iteration > maxNumOfIteration){
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

//    DIISprocedure();

    eig_sym(eigVal, eigVec, V.t()*m_Fu*V);
    m_Cu = V*eigVec;
    m_fockEnergyU = eigVal;

    eig_sym(eigVal, eigVec, V.t()*m_Fd*V);
    m_Cd = V*eigVec;
    m_fockEnergyD = eigVal;


    normalize(m_Cu, m_nSpinUpElectrons);
    normalize(m_Cd, m_nSpinDownElectrons);

    m_Pu = 0.5 * m_Pu + 0.5*m_Cu.cols(0, m_nSpinUpElectrons-1) * m_Cu.cols(0, m_nSpinUpElectrons-1).t();
    m_Pd = 0.5 * m_Pd + 0.5*m_Cd.cols(0, m_nSpinDownElectrons-1) * m_Cd.cols(0, m_nSpinDownElectrons-1).t();
}


void UHF::DIISprocedure()
{
    m_errorsU.push_back(m_Fu*m_Pu*m_S - m_S*m_Pu*m_Fu);
    m_fockMatricesU.push_back(m_Fu);

    m_errorsD.push_back(m_Fd*m_Pd*m_S - m_S*m_Pd*m_Fd);
    m_fockMatricesD.push_back(m_Fd);


    if(m_errorsU.size() > 20){
        m_errorsU.erase(m_errorsU.begin());
        m_fockMatricesU.erase(m_fockMatricesU.begin());

        m_errorsD.erase(m_errorsD.begin());
        m_fockMatricesD.erase(m_fockMatricesD.begin());

        mat Au = zeros(m_errorsU.size()+1, m_errorsU.size()+1);
        mat Ad = zeros(m_errorsD.size()+1, m_errorsD.size()+1);

        for(uint i = 0; i < m_errorsU.size(); i++){
            Au(i, m_errorsU.size()) = -1;
            Au(m_errorsU.size(), i) = -1;

            Ad(i, m_errorsD.size()) = -1;
            Ad(m_errorsD.size(), i) = -1;

            for(uint j = 0; j < m_errorsU.size(); j++){
                Au(i,j) = trace(m_errorsU.at(i) * m_errorsU.at(j));
                Ad(i,j) = trace(m_errorsD.at(i) * m_errorsD.at(j));
            }
        }

        vec bu = zeros(Au.n_rows);
        vec bd = zeros(Ad.n_rows);

        bu(bu.n_elem -1) = -1.0;
        bd(bd.n_elem -1) = -1.0;
        bu = inv(Au) * bu;
        bd = inv(Ad) * bd;

        for(uint i = 0; i < bu.n_elem - 1; i++){
            m_Fu += bu(i) * m_fockMatricesU.at(i);
            m_Fd += bd(i) * m_fockMatricesD.at(i);
        }
    }
}



void UHF::calculateEnergy()
{
    m_energy = 0.5 * accu( (m_Pu + m_Pd) % m_h + m_Fu % m_Pu + m_Fd % m_Pd)
               + m_system->getNucleiPotential();
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

field<mat> UHF::getFockMatrix()
{
    updateFockMatrix();
    field<mat> fockMatrices(2,1);
    fockMatrices(0) = m_Fu;
    fockMatrices(1) = m_Fd;
    return fockMatrices;
}

field<mat> UHF::getDensityMatrix() const
{
    field<mat> densityMatrices(2,1);
    densityMatrices(0) = m_Pu;
    densityMatrices(1) = m_Pd;
    return densityMatrices;
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
                    sumDensity += (m_Pu(p,p) + m_Pd(p,p)) * innerProduct * dr;
                    m_density(j,i,k) += (m_Pu(p,p) + m_Pd(p,p)) * innerProduct ;

                    for(int q = p+1; q < m_nBasisFunctions; q++){
                        innerProduct = m_system->gaussianProduct(p, q, x(i), y(j), z(k));
                        sumDensity += 2.0 * (m_Pu(p,p) + m_Pd(p,p)) * innerProduct * dr;
                        m_density(j,i,k) += 2.0 * (m_Pu(p,p) + m_Pd(p,p)) * innerProduct ;

                    }
                }

            }
        }
    }

    cout << "density sum: " << sumDensity << endl;
    densityOutput(x.min(),x.max(),y.min(),y.max(),z.min(),z.max());

}
