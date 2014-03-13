#include "bomd.h"

using namespace hf;

BOMD::BOMD(System *system, HFsolver *solver, const int &rank, const int &nProcs):
    m_rank(rank),
    m_nProcs(nProcs),
    m_system(system),
    m_solver(solver),
    m_nCores(system->getNumOfCores()),
    m_nOrbitals(system->getTotalNumOfBasisFunc())


{
    m_nSteps = 100;
    m_dtn   =  4.0;
    m_dampingFactor = 0.0;


    pos    = zeros(m_nCores, 3);
    posNew = zeros(m_nCores, 3);
    posOld = zeros(m_nCores, 3);

    for(int core = 0; core < m_nCores; core++){
        pos.row(core) = m_system->m_basisSet.at(core)->corePosition();
    }

    posOld = pos;
    m_GD = new GeometricalDerivative(m_system,m_solver);
}


void BOMD::runDynamics()
{
    for(int nStep = 0; nStep < m_nSteps; nStep++){
        solveSingleStep();
        writeToFile(pos,nStep);
        posOld = pos;
        pos = posNew;
        updateCorePositions();

    }

}

void BOMD::solveSingleStep()
{
    m_solver->runSolver();
    m_energy  = m_solver->getEnergy();


    for(int core = 0; core < m_nCores; core++){
        clock_t begin = clock();
        m_energyGradient = m_GD->energyGradient(core);
        clock_t end = clock();

        if(m_rank==0){
            cout << setprecision(3)
                 << "Elapsed time on gradient calculation: "
                 << (double(end - begin))/CLOCKS_PER_SEC << endl;
        }

        IntegrateCoreForwardInTime(core);
    }


    if(m_rank==0){
        cout << "-------------------------------------------------------------------------------------"  << endl;
    }
}



void BOMD::IntegrateCoreForwardInTime(int core)
{
    int coreMass = m_system->m_basisSet.at(core)->coreMass();

    posNew.row(core) = 2 * pos.row(core) - posOld.row(core)
            - m_dtn * m_dtn * m_energyGradient / coreMass;


}



void BOMD::updateCorePositions()
{
    for(int core = 0; core < m_nCores; core++){
        m_system->m_basisSet.at(core)->setCorePosition(pos.row(core));
    }

}

double BOMD::getEnergy() const
{
    return m_energy;
}
const rowvec& BOMD::getEnergyGradient() const
{
    return m_energyGradient;
}




void BOMD::writeToFile(mat R, int currentTimeStep) {
    ivec atomTypes;

    if(m_nCores == 3){
        atomTypes << 1 << 1 << 8 << 1;
    }else if(m_nCores == 5){
         atomTypes << 14 << 8 << 8 << 8 << 8;
    }
    else{
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

