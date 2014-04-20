#include "bomd.h"

using namespace hf;

BOMD::BOMD(ElectronicSystem *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_atoms(system->atoms()),
    m_nAtoms(system->nAtoms())
{
    m_nSteps = 500;
    m_dt   =  0.1;
    m_frictionConstant = 0.4;

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
    computeForces();
    for(int i = 0; i < m_nSteps; i++){
        if(m_rank == 0 ){
            cout << "MD step:   " << i << endl;
            cout << "-------------------------------------------------------------------------------------"  << endl;

        }

        solveSingleStep();
        writeLammpsFile(i);
        systemProperties(i);
    }
    writeSystemProperties();
}


void BOMD::computeForces()
{
    m_solver->runSolver();
    m_energyGradient = -m_GD->energyGradient();
}

void BOMD::solveSingleStep()
{
    halfKick();
    for(Atom* atom : m_atoms){
        rowvec corePosition = atom->corePosition() + m_dt * atom->coreVelocity();
        atom->setCorePosition(corePosition);
    }
    computeForces();
    halfKick();
}

void BOMD::halfKick()
{
    int i = 0;
    for(Atom* atom : m_atoms){
        rowvec coreVelocity = atom->coreVelocity() + 0.5 * m_dt
                            * (m_energyGradient.row(i)/atom->coreMass()
                            -  m_frictionConstant * atom->coreVelocity());

        atom->setCoreVelocity(coreVelocity);
        i++;
    }

}


void BOMD::systemProperties(int currentTimeStep)
{
    int i = currentTimeStep;
    m_time(i)         = i * m_dt;
    m_potentialEnergy(i) = m_solver->energy();
    for (const Atom* atom : m_atoms){
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

    stringstream outStepName;
    outStepName <<"/home/milad/kurs/qmd/state" << setw(4) << setfill('0')  << currentTimeStep <<".lmp";
    ofstream lammpsFile(outStepName.str(), ios::out | ios::binary);

    // The system boundaries
    double xMin = -5.0;
    double xMax = 5.0;
    double yMin = -5.0;
    double yMax = 5.0;
    double zMin = -5.0;
    double zMax = 5.0;
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
        double atomType = m_atoms[i]->atomType();
        rowvec position = m_atoms[i]->corePosition();
        rowvec velocity = m_atoms[i]->coreVelocity();

        lammpsFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&position(2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&velocity(2)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,0)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,1)), sizeof(double));
        lammpsFile.write(reinterpret_cast<const char*>(&m_energyGradient(i,2)), sizeof(double));

    }
    lammpsFile.close();
}

