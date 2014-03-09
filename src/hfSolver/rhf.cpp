#include "rhf.h"

using namespace hf;

RHF::RHF(System *system, const int &rank, const int &nProcs):
    HFsolver(system, rank, nProcs),
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
    int step = 0;
    int maxStep = 1000;

    while (stdDeviation  > HFSOLVERTOLERANCE){
        fockEnergyOld = m_fockEnergy;
        solveSingle();
        stdDeviation = computeStdDeviation(m_fockEnergy, fockEnergyOld);
        updateFockMatrix();

        step+=1;
        if(step > maxStep){
            cerr << "Energy has not converged! " << endl;
            exit(1);
        }
    }

}
void RHF::solveSingle()
{
    vec eigVal;
    mat eigVec;
    eig_sym(eigVal, eigVec, m_S);

    mat V = eigVec*diagmat(1.0/sqrt(eigVal));

    eig_sym(eigVal, eigVec, V.t() * m_F * V);

    m_C = V*eigVec;

    normalize(m_C, m_nElectrons/2);

    m_P =  2.0 * m_C.cols(0, m_nElectrons/2.0-1) * m_C.cols(0, m_nElectrons/2.0-1).t();
    m_fockEnergy = eigVal;
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
    m_energy = 0.5 * accu(m_P % (m_h + m_F)) + m_system->getNucleiPotential();
}

const mat& RHF::getExpansionCoeff() const
{
    return m_C;
}


field<mat> RHF::getFockMatrix()
{
    updateFockMatrix();
    field<mat> fockMatrices(1,1);
    fockMatrices(0,0) = m_F;
    return fockMatrices;
}

field<mat> RHF::getDensityMatrix() const
{
    field<mat> densityMatrices(1,1);
    densityMatrices(0,0) = m_P;
    return densityMatrices;
}
void RHF::calculateDensity()
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
                    sumDensity += m_P(p,p) * innerProduct * dr;
                    m_density(j,i,k) += m_P(p,p) * innerProduct ;

                    for(int q = p+1; q < m_nBasisFunctions; q++){
                        innerProduct = m_system->gaussianProduct(p, q, x(i), y(j), z(k));
                        sumDensity += 2.0 * m_P(p,q) * innerProduct * dr;
                        m_density(j,i,k) += 2.0 * m_P(p,q) * innerProduct ;

                    }
                }

            }
        }
    }


    //    for(int i = xMin; i < xMax; i++) {
    //        for(int j = yMin; j < yMax; j++) {
    //            for(int k = zMin; k < zMax; k++) {

    //                for(int p = 0; p < m_nBasisFunctions; p++){
    //                    for(int q = 0; q < m_nBasisFunctions; q++){
    //                        double innerProduct = m_system->gaussianProduct(p, q, x(i), y(j), z(k));
    //                        m_density(j,i,k) += 2.0 * m_C(p,6)* m_C(q,6) * innerProduct ;

    //                    }
    //                }

    //            }
    //        }
    //    }




    cout << "density sum: " << sumDensity << endl;
    //    cout << m_density << endl;
    //    cout <<"rank: "<<m_rank<<" xLim: " <<xMin << "   "<< xMax << endl;
    //    cout <<"rank: "<<m_rank<<" yLim: "<< yMin << "   "<< yMax << endl;
    //    cout <<"rank: "<<m_rank<<" zLim: "<< zMin << "   "<< zMax << endl;
    densityOutput(x.min(),x.max(),y.min(),y.max(),z.min(),z.max());
}
