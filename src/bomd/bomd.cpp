#include "bomd.h"

using namespace hf;

BOMD::BOMD(ElectronicSystem *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_atoms(system->atoms()),
    m_nAtoms(system->nAtoms()),
    m_corePositions(zeros(m_nAtoms, 3)),
    m_coreVelocities(zeros(m_nAtoms, 3))
{
    m_nSteps = 20;
    m_dt   =  0.1;
    m_dampingFactor = 0.0;

    m_time = zeros(m_nSteps);
    m_kineticEnergy = zeros(m_nSteps);
    m_potentialEnergy = zeros(m_nSteps);
    m_totalEnergy = zeros(m_nSteps);

    m_GD = new GeometricalDerivative(m_system, m_solver);

    m_rank = 0;
#if USE_MPI
    boost::mpi::communicator world;
    m_rank = world.rank();
#endif

}
void BOMD::runDynamics()
{
    initialStep();
    for(int i = 0; i < m_nSteps; i++){
        if(m_rank == 0 ){
            cout << "MD step:   " << i << endl;
            cout << "-------------------------------------------------------------------------------------"  << endl;

        }

        solveSingleStep();
        writeLammpsFile(i);
        updateCores();
        systemProperties(i);
    }
    writeSystemProperties();
}

void BOMD::initialStep()
{
    for(int i = 0; i < m_nAtoms; i++){
        m_corePositions.row(i)  = m_atoms[i]->corePosition();
        m_coreVelocities.row(i) = m_atoms[i]->coreVelocity();
    }
    computeForces();
}


void BOMD::computeForces()
{
    m_solver->runSolver();
    m_energyGradient = -m_GD->energyGradient();
}


void BOMD::updateCores()
{
    for(int i = 0; i < m_nAtoms; i++){
        m_atoms[i]->setCorePosition(m_corePositions.row(i));
        m_atoms[i]->setCoreVelocity(m_coreVelocities.row(i));
    }
}

void BOMD::solveSingleStep()
{
    halfKick();
    for (int i=0; i< m_nAtoms; i++){
        m_corePositions.row(i) += m_dt * m_coreVelocities.row(i);
    }
    updateCores();
    computeForces();
    halfKick();
}

void BOMD::halfKick()
{
    for (int i=0; i < m_nAtoms; i++){
        m_coreVelocities.row(i) += 0.5 * m_dt * m_energyGradient.row(i)/m_atoms[i]->coreMass();
    }
}


void BOMD::systemProperties(int currentTimeStep)
{
    int i = currentTimeStep;
    m_time(i)         = i * m_dt;
    m_potentialEnergy(i) = m_solver->energy();
    for (Atom* atom : m_atoms){
        m_kineticEnergy(i) += 0.5 * atom->coreMass() * dot(atom->coreVelocity(),atom->coreVelocity());
    }

    m_totalEnergy(i) = m_kineticEnergy(i) + m_potentialEnergy(i);
}


double BOMD::potentialEnergy() const
{
    return m_solver->energy();
}

const mat& BOMD::energyGradient() const
{
    return m_energyGradient;
}


void BOMD::writeSystemProperties()
{
    fstream output;
    output.open ( "/home/milad/kurs/qmd/systemProperties/properties.txt",fstream::out);
    output  << "Time    "  <<"Kinetic       "
            << "Potential     " <<"Total Energy     "
            <<endl;

    for(int state = 0; state < m_nSteps; state++){
        output << m_time(state)
               <<"       "<< m_kineticEnergy(state)
               <<"       "<< m_potentialEnergy(state)
               <<"       "<< m_totalEnergy(state)
               << endl;
    }
    output.close();
}


void BOMD::writeLammpsFile(int currentTimeStep) {
    ivec atomTypes;

    if(m_nAtoms == 3){
        atomTypes << 8 << 1 << 1;
    }else if(m_nAtoms == 5){
         atomTypes << 14 << 8 << 8 << 8 << 8;
    }else if(m_nAtoms == 9){
        atomTypes << 14 << 8 << 8 << 8 << 8 << 14 << 8 << 8 <<8;
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
    int nColumns = 1 + 3 + 3 + 3;
    // We could divide the data into chunks by the LAMMPS file format, but we don't - i.e. only
    // use one chunk. The chunk length is then the same as the number of atoms times the number
    // of columns.
    int nChunks = 1;
    int chunkLength = m_nAtoms * nColumns;

    // Write all the above to the lammps file
    lammpsFile.write(reinterpret_cast<const char*>(&currentTimeStep), sizeof(int));
    lammpsFile.write(reinterpret_cast<const char*>(&m_nAtoms), sizeof(int));
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
    for(int i = 0; i < m_nAtoms; i++) {
        // IMPORTANT: Even though atom numbers are usually integers, they must be written
        // as double according to the LAMMPS standard.
        double atomType = atomTypes(i);
        lammpsFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));

        // Write the x, y and z-components
        lammpsFile.write(reinterpret_cast<const char*>(&m_corePositions(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_corePositions(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_corePositions(i,2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_coreVelocities(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_coreVelocities(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_coreVelocities(i,2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,2)), sizeof(double));

    }
    lammpsFile.close();
}

