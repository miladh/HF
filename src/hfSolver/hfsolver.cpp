#include "hfsolver.h"

using namespace hf;

HFsolver::HFsolver(System *system, const int &rank, const int &nProcs):
    m_rank(rank),
    m_nProcs(nProcs),
    m_step(0),
    m_system(system),

    m_nBasisFunctions(system->getTotalNumOfBasisFunc()),
    m_S(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_h(zeros(m_nBasisFunctions,m_nBasisFunctions))

{
    m_Q.set_size(m_nBasisFunctions, m_nBasisFunctions);
    for(int i = 0; i < m_nBasisFunctions; i++){
        for(int j = 0; j < m_nBasisFunctions; j++){
            m_Q(i,j) = zeros(m_nBasisFunctions,m_nBasisFunctions);
        }
    }

}


void HFsolver::runSolver()
{

    setupOneParticleMatrix();
    setupTwoParticleMatrix();
    updateFockMatrix();
    advance();
    calculateEnergy();

    cout << setprecision(14)
         << "configuration " << m_step
         << " - Energy: "  << m_energy << endl;

    //    calculateDensity();
    m_step+=1;
}


void HFsolver::setupOneParticleMatrix()
{
    rowvec oneElectronIntegrals;

    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int q = p; q < m_nBasisFunctions; q++){
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
    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int r = 0; r < m_nBasisFunctions; r++){
            for(int q = p; q < m_nBasisFunctions; q++){
                for(int s = r; s < m_nBasisFunctions; s++){

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



const double& HFsolver::getEnergy() const
{
    return m_energy;
}

const field<mat>& HFsolver::getQmatrix(){

    return m_Q;
}
const mat& HFsolver::gethmatrix(){

    return m_h;
}

const mat& HFsolver::getSmatrix(){

    return m_S;
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


