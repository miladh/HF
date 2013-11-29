#include "bomd.h"

BOMD::BOMD(System *system):
    m_system(system),
    m_nCores(system->getNumOfCores()),
    m_nElectrons(system->getNumOfElectrons()),
    m_nOrbitals(system->getTotalNumOfBasisFunc())


{
    m_nSteps = 300;
    m_dtn   =  4.0;
    m_dampingFactor = 0.0;

    //initialize:
    m_dh.set_size(m_nOrbitals,m_nOrbitals);
    m_dS.set_size(m_nOrbitals,m_nOrbitals);
    m_dQ.set_size(m_nOrbitals, m_nOrbitals);
    for(int i = 0; i < m_nOrbitals; i++ ){
        for(int j = 0; j < m_nOrbitals; j++ ){
            m_dh(i,j)  = zeros<rowvec>(3);
            m_dS(i,j)  = zeros<rowvec>(3);
            m_dQ(i,j).set_size(m_nOrbitals,m_nOrbitals);
        }
    }

    for(int i = 0; i < m_nOrbitals; i++ ){
        for(int j = 0; j < m_nOrbitals; j++ ){
            for(int k = 0; k < m_nOrbitals; k++ ){
                for(int l = 0; l < m_nOrbitals; l++ ){
                    m_dQ(i,j)(k,l) = zeros<rowvec>(3);
                }
            }
        }
    }

    pos    = zeros(m_nCores, 3);
    posNew = zeros(m_nCores, 3);
    posOld = zeros(m_nCores, 3);

    for(int core = 0; core < m_nCores; core++){
        pos.row(core) = m_system->m_basisSet.at(core)->corePosition();
    }

    posOld = pos;
    m_solver = new HFsolver(m_system);
}


void BOMD::runDynamics()
{

    for(int nStep = 0; nStep < m_nSteps; nStep++){
        solveSingleStep();
        writeToFile(pos,nStep);
        posOld = pos;
        pos = posNew;
        updateCorePositions();

        cout << "step " << nStep << " Energy: " << m_energy << endl;
    }

}

void BOMD::solveSingleStep()
{
    m_solver->runSolver();
    m_P = m_solver->getDensityMatrix();
    m_energy  = m_solver->getEnergy();

    for(int core = 0; core < m_nCores; core++){
        setupDerivativeMatrices(core);
        m_energyGradient = calculateEnergyGradient(core);
        IntegrateCoreForwardInTime(core);
    }
}




void BOMD::IntegrateCoreForwardInTime(int core)
{
    int coreMass = m_system->m_basisSet.at(core)->coreMass();

    posNew.row(core) = 2 * pos.row(core) - posOld.row(core)
            - m_dtn * m_dtn * m_energyGradient / coreMass;


}

rowvec BOMD::calculateEnergyGradient(int core)
{
    rowvec dE  = {0,0,0};

    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){
            dE += m_P(p, q)*m_dh(p, q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    dE += 0.5*m_P(p,q)*m_P(s,r)*(m_dQ(p,r)(q,s) - 0.5*m_dQ(p,r)(s,q));
                }
            }
        }
    }


    mat dSx, dSy,dSz;
    dSx = zeros(m_nOrbitals,m_nOrbitals);
    dSy = zeros(m_nOrbitals,m_nOrbitals);
    dSz = zeros(m_nOrbitals,m_nOrbitals);
    for(int i = 0; i < m_nOrbitals; i++){
        for(int j = 0; j < m_nOrbitals; j++){
            dSx(i,j) = m_dS(i, j)(0);
            dSy(i,j) = m_dS(i, j)(1);
            dSz(i,j) = m_dS(i, j)(2);
        }
    }
    mat F = m_solver->getF();

    dE(0) -= 0.5 * trace(m_P*dSx*m_P*F);
    dE(1) -= 0.5 * trace(m_P*dSy*m_P*F);
    dE(2) -= 0.5 * trace(m_P*dSz*m_P*F);

    //    Nuclear repulsion term
    dE  +=m_system->getNucleiPotential_derivative(core);

    return dE;
}


void BOMD::setupDerivativeMatrices(const int core)
{
    mat diffOneParticleIntegral;
    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            diffOneParticleIntegral = m_system->getOneParticleDerivative(p,q,core);
            m_dS(p,q) = diffOneParticleIntegral.row(0);
            m_dh(p,q) = diffOneParticleIntegral.row(1);
        }
    }


    //Set up the dQ array:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int r = 0; r < m_nOrbitals; r++){
            for(int q = 0; q < m_nOrbitals; q++){
                for(int s = 0; s < m_nOrbitals; s++){

                    m_dQ(p,r)(q,s) = m_system->getTwoParticleIntegralDerivative(p,q,r,s,core);

                }
            }
        }
    }

}


void BOMD::updateCorePositions()
{

    for(int core = 0; core < m_nCores; core++){
        m_system->m_basisSet.at(core)->setCorePosition(pos.row(core));
    }

}

rowvec BOMD::getEnergyGradient() const
{
    return m_energyGradient;
}

double BOMD::getEnergy() const
{
    return m_energy;
}


void BOMD::writeToFile(mat R, int currentTimeStep) {
    ivec atomTypes;

    if(m_nCores > 2){
        atomTypes << 1 << 1 << 8 << 1;
    }else{
        atomTypes << 1 << 1;
    }

    stringstream outStepName;
    outStepName <<"/home/milad/kurs/qmd/state" << setw(4) << setfill('0')  << currentTimeStep <<".lmp";
    ofstream lammpsFile(outStepName.str(), ios::out | ios::binary);


    // The system boundaries
    double xMin = -2.0;
    double xMax = 2.0;
    double yMin = -2.0;
    double yMax = 2.0;
    double zMin = -2.0;
    double zMax = 2.0;
    // Shearing is zero unless the system boundaries are sheared (yes that's "sheared",
    // not "shared")
    double xShear = 0.0;
    double yShear = 0.0;
    double zShear = 0.0;
    // nColumns is the number of data types you want to write. In our case we want to
    // write four - the atom type and the x, y and z components of the position.
    // If you want velocities, forces, etc., just add more columns and write more data.
    int nColumns = 1 + 3 + 2;
    // We could divide the data into chunks by the LAMMPS file format, but we don't - i.e. only
    // use one chunk. The chunk length is then the same as the number of atoms times the number
    // of columns.
    int nChunks = 1;
    int chunkLength = R.n_rows * nColumns;

    // Write all the above to the lammps file
    lammpsFile.write(reinterpret_cast<const char*>(&currentTimeStep), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&R.n_rows), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&xMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&xMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zMin), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zMax), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&xShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&yShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&zShear), sizeof(double));
    lammpsFile.write(reinterpret_cast<const char*>(&nColumns), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&nChunks), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&chunkLength), sizeof(int));

    // Write all the data for each atom to file
    for(uint i = 0; i < R.n_rows; i++) {
        // IMPORTANT: Even though atom numbers are usually integers, they must be written
        // as double according to the LAMMPS standard.
        double atomType = atomTypes(i);
        double force = m_energyGradient(0);
        lammpsFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));

        // Write the x, y and z-components
        lammpsFile.write(reinterpret_cast<const char*>(&R(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&R(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&R(i,2)), sizeof(double));

        lammpsFile.write(reinterpret_cast<const char*>(&m_energy), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&force), sizeof(double));

    }
    lammpsFile.close();
}
