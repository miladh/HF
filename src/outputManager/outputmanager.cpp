#include "outputmanager.h"


using namespace hf;
using namespace H5;


OutputManager::OutputManager():
    m_rank(0),
    m_nProcs(0)
{



#if USE_MPI
    boost::mpi::environment env;
    boost::mpi::communicator world;
    m_rank = world.rank();
    m_nProcs = world.size();
#endif

    m_outputFileName << "/home/milad/kurs/qmd/output_" << m_rank << ".h5";
    m_output = new H5File (m_outputFileName.str(), H5F_ACC_TRUNC);
}

void OutputManager::registerAtoms(vector<Atom *> atoms)
{

    CompType atomCompound( sizeof(AtomProperties) );
     atomCompound.insertMember( "type", HOFFSET(AtomProperties, type), PredType::NATIVE_INT);
     atomCompound.insertMember( "basisType", HOFFSET(AtomProperties, basisType), StrType(PredType::C_S1, 64));
     atomCompound.insertMember("x", HOFFSET(AtomProperties, x), PredType::NATIVE_DOUBLE);
     atomCompound.insertMember("y", HOFFSET(AtomProperties, y), PredType::NATIVE_DOUBLE);
     atomCompound.insertMember("z", HOFFSET(AtomProperties, z), PredType::NATIVE_DOUBLE);
     atomCompound.insertMember("partialCharge", HOFFSET(AtomProperties, corePartialCharge), PredType::NATIVE_DOUBLE);
     atomCompound.insertMember("coreCharge", HOFFSET(AtomProperties, coreCharge), PredType::NATIVE_DOUBLE);

    hsize_t dim[1];
    dim[0] = atoms.size();

    DataSpace space(1, dim);
    DataSet *dataset = new DataSet(m_output->createDataSet("DatasetName", atomCompound, space));

    AtomProperties * atomProp = new AtomProperties[atoms.size()];



    for(int i = 0; i < signed(atoms.size()); i++) {
        Atom* atom = atoms.at(i);
        strcpy(atomProp[i].basisType, atom->basisType().c_str());
        atomProp[i].type = atom->atomType();
        atomProp[i].x = 4;
        atomProp[i].y = 5;
        atomProp[i].z = 7;
        atomProp[i].coreCharge = atom->coreCharge();
        atomProp[i].corePartialCharge = atom->corePartialCharge();
        dataset->write(atomProp, atomCompound);
    }

    double energy = 2;
    Attribute energyAttribute(dataset->createAttribute("energy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR));
    energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);


    m_output->flush(H5F_SCOPE_GLOBAL);
    m_output->close();

}

void OutputManager::saveEnergy(const double& energy)
{
//    Attribute energyAttribute(stateDataSet.createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
//    energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);

}

