#include "hfsolver.h"

HFsolver::HFsolver(System *system, const int &rank, const int &nProcs):
    m_rank(rank),
    m_nProcs(nProcs),
    m_system(system),
    m_nElectrons(system->getNumOfElectrons()),
    m_nOrbitals(system->getTotalNumOfBasisFunc()),
    m_S(zeros(m_nOrbitals,m_nOrbitals)),
    m_h(zeros(m_nOrbitals,m_nOrbitals)),
    m_F(zeros(m_nOrbitals,m_nOrbitals)),
    m_P(zeros(m_nOrbitals,m_nOrbitals)),
    m_C(ones(m_nOrbitals,m_nOrbitals)),
    m_step(0)

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
//    calculateDensity();
    m_step+=1;
}

void HFsolver::solveSingle()
{
    vec eigVal;
    mat eigVec;
    eig_sym(eigVal, eigVec, m_S);

    mat V = eigVec*diagmat(1.0/sqrt(eigVal));

    m_F = V.t()*m_F*V;


    eig_sym(eigVal, eigVec, m_F);
    m_C = V*eigVec;

    normalize();


//    m_P = 0.5*m_P + m_C * m_C.t();  // Interpolate between new and old density matrix
    m_P = 2*m_C.cols(0, m_nElectrons/2.0-1)*m_C.cols(0, m_nElectrons/2.0-1).t();

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

    cout << setprecision(14) << "configuration " << m_step << " Energy: "  << m_energy << endl;



}

void HFsolver::setupOneParticleMatrix()
{
    rowvec oneElectronIntegrals;

    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = p; q < m_nOrbitals; q++){
            oneElectronIntegrals = m_system->getOneParticleIntegral(p,q);
            m_S(p,q) = oneElectronIntegrals(0);
            m_h(p,q) = oneElectronIntegrals(1);
        }
    }

    m_S = symmatu(m_S);
    m_h = symmatu(m_h);


}


void HFsolver::setupTwoParticleMatrix()
{
    for(int p = 0; p < m_nOrbitals; p++){
        for(int r = 0; r < m_nOrbitals; r++){
            for(int q = p; q < m_nOrbitals; q++){
                for(int s = r; s < m_nOrbitals; s++){

                    m_Q(p,r)(q,s) = m_system->getTwoParticleIntegral(p,q,r,s);
                    m_Q(q,r)(p,s) = m_Q(p,r)(q,s);
                    m_Q(p,s)(q,r) = m_Q(p,r)(q,s);
                    m_Q(q,s)(p,r) = m_Q(p,r)(q,s);
                    m_Q(r,p)(s,q) = m_Q(p,r)(q,s);
                    m_Q(s,p)(r,q) = m_Q(p,r)(q,s);
                    m_Q(r,q)(s,p) = m_Q(p,r)(q,s);
                    m_Q(s,q)(r,p) = m_Q(p,r)(q,s);
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


void HFsolver::calculateDensity()
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

                for(int p = 0; p < m_nOrbitals; p++){
                    double innerProduct = m_system->gaussianProduct(p, p, x(i), y(j), z(k));
                    sumDensity += m_P(p,p) * innerProduct * dr;
                    m_density(j,i,k) += m_P(p,p) * innerProduct ;

                    for(int q = p+1; q < m_nOrbitals; q++){
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

//                for(int p = 0; p < m_nOrbitals; p++){
//                    for(int q = 0; q < m_nOrbitals; q++){
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


void HFsolver::densityOutput(const double &xMin, const double &xMax,
                             const double &yMin, const double &yMax,
                             const double &zMin, const double &zMax)
{

    stringstream cubeFileName;
    cubeFileName<<"/home/milad/kurs/qmd/density/id"<< m_rank <<"_cubeFile" << setw(4) << setfill('0')  << m_step <<".bin";
    ofstream cubeFile(cubeFileName.str(), ios::out | ios::binary);

    //Header
    double nCores = m_system->getNumOfCores();
    vec origo = {0,0,0};
    vec xLimits = {xMin, xMax};
    vec yLimits = {yMin, yMax};
    vec zLimits = {zMin, zMax};
    double nX = m_density.n_rows;
    double nY = m_density.n_cols;
    double nZ = m_density.n_slices;

    cubeFile.write(reinterpret_cast<const char*>(&nCores), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(1)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(2)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nX), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&xLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&xLimits(1)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nY), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&yLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&yLimits(1)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nZ), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&zLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&zLimits(1)), sizeof(double));

    for(int core = 0; core < nCores; core++){
        rowvec3 R = m_system->m_basisSet.at(core)->corePosition();
        double atomType = m_system->m_basisSet.at(core)->coreCharge();

        cubeFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));

        // Write the x, y and z-components
        cubeFile.write(reinterpret_cast<const char*>(&R(0)), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&R(1)), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&R(2)), sizeof(double));
    }

    //density data
    for(uint k = 0; k < m_density.n_slices; k++) {
        for(uint i = 0; i < m_density.n_rows; i++) {
            for(uint j = 0; j < m_density.n_cols; j++) {

                cubeFile.write(reinterpret_cast<const char*>(&m_density(i,j,k)), sizeof(double));
            }
        }
    }

    cubeFile.close();
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

double HFsolver::getFockEnergy() const
{
     return m_fockEnergy;
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

mat HFsolver::getF()
{
    setupFockMatrix();
    return m_F;
}

mat HFsolver::getDensityMatrix() const
{
    return m_P;
}







