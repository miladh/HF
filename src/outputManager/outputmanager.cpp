#include "outputmanager.h"


using namespace hf;

OutputManager::OutputManager(const Config* cfg, const int nAtoms):
    m_cfg(cfg),
    m_rank(0),
    m_nProcs(0)
{

#if USE_MPI
    boost::mpi::environment env;
    boost::mpi::communicator world;
    m_rank = world.rank();
    m_nProcs = world.size();
#endif

    const Setting & root = m_cfg->getRoot();
    string outputFilePath = root["analysisSettings"]["outputFilePath"];
    m_outputFileName << outputFilePath << "/HFOutput_" << m_rank << ".h5";
    m_output = new H5File (m_outputFileName.str(), H5F_ACC_TRUNC);


    Group rootGroup = m_output->openGroup("/");
    Attribute nAtoms_a(rootGroup.createAttribute("nAtoms", PredType::NATIVE_INT, H5S_SCALAR));
    nAtoms_a.write(PredType::NATIVE_INT, &nAtoms);

    //---------------------------------------------------------------------------------------------------------
    m_atomCompound = new CompType(sizeof(AtomAttributes));
    m_atomCompound->insertMember( "type", HOFFSET(AtomAttributes, type), PredType::NATIVE_INT);
    m_atomCompound->insertMember( "basis type", HOFFSET(AtomAttributes, basisType), StrType(PredType::C_S1, 64));
    m_atomCompound->insertMember("x", HOFFSET(AtomAttributes, x), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("y", HOFFSET(AtomAttributes, y), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("z", HOFFSET(AtomAttributes, z), PredType::NATIVE_DOUBLE);
    m_atomCompound->insertMember("core charge", HOFFSET(AtomAttributes, coreCharge), PredType::NATIVE_INT);
    m_atomCompound->insertMember("partial charge", HOFFSET(AtomAttributes, corePartialCharge), PredType::NATIVE_DOUBLE);

   hsize_t dim[1];
   dim[0] = nAtoms;
   DataSpace space(1, dim);
   m_dataset = new DataSet(m_output->createDataSet("atoms", *m_atomCompound, space));
   m_atomAttributes = new AtomAttributes[nAtoms];
   //---------------------------------------------------------------------------------------------------------
}

void OutputManager::saveAtoms(vector<Atom *> atoms)
{

    for(int i = 0; i < signed(atoms.size()); i++) {
        Atom* atom = atoms.at(i);
        cout << atom->basisType() << endl;
        strcpy(m_atomAttributes[i].basisType, atom->basisType().c_str());
        m_atomAttributes[i].type = atom->atomType();
        m_atomAttributes[i].x = atom->corePosition()(0);
        m_atomAttributes[i].y = atom->corePosition()(1);
        m_atomAttributes[i].z = atom->corePosition()(2);
        m_atomAttributes[i].coreCharge = atom->coreCharge();
        m_atomAttributes[i].corePartialCharge = atom->corePartialCharge();
        m_dataset->write(m_atomAttributes, *m_atomCompound);
    }


}

void OutputManager::saveEnergy(const double& energy, const mat& orbitalEnergies)
{

    Attribute energyAttribute(m_dataset->createAttribute("energy", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    energyAttribute.write(PredType::NATIVE_DOUBLE, &energy);

    hsize_t dim[2] = {orbitalEnergies.n_cols, orbitalEnergies.n_rows};
    DataSpace space(2, dim);
    DataSet dataset(m_output->createDataSet("orbital energies", PredType::NATIVE_DOUBLE, space));
    dataset.write(orbitalEnergies.memptr(), PredType::NATIVE_DOUBLE);

    if(orbitalEnergies.n_rows > 1){
        string comment = " Col 0 is spin up, col 1 is spin down";
        StrType vlst(0, H5T_VARIABLE);
        Attribute attrConfigFilename = dataset.createAttribute("comment", vlst, DataSpace(H5S_SCALAR));
        attrConfigFilename.write(vlst, comment);
    }
}

void OutputManager::saveElectronDensity(const field<cube>& densityCubes)
{    
    const field<cube>& density = densityCubes;
    Group electronDensities( m_output->createGroup( "/electronDensity" ));

    const Setting & root = m_cfg->getRoot();
    int    nGridPoints  = root["analysisSettings"]["densitySettings"]["nGridPoints"];
    double minGridPoint = root["analysisSettings"]["densitySettings"]["minGridPoint"];
    double maxGridPoint = root["analysisSettings"]["densitySettings"]["maxGridPoint"];

    Attribute nXPoints(electronDensities.createAttribute("nXPoints", PredType::NATIVE_INT, H5S_SCALAR));
    Attribute xMin(electronDensities.createAttribute("xMin", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute xMax(electronDensities.createAttribute("xMax", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    nXPoints.write(PredType::NATIVE_INT, &nGridPoints);
    xMin.write(PredType::NATIVE_DOUBLE, &minGridPoint);
    xMax.write(PredType::NATIVE_DOUBLE, &maxGridPoint);

    Attribute nYPoints(electronDensities.createAttribute("nYPoints", PredType::NATIVE_INT, H5S_SCALAR));
    Attribute yMin(electronDensities.createAttribute("yMin", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute yMax(electronDensities.createAttribute("yMax", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    nYPoints.write(PredType::NATIVE_INT, &nGridPoints);
    yMin.write(PredType::NATIVE_DOUBLE, &minGridPoint);
    yMax.write(PredType::NATIVE_DOUBLE, &maxGridPoint);

    Attribute nZPoints(electronDensities.createAttribute("nZPoints", PredType::NATIVE_INT, H5S_SCALAR));
    Attribute zMin(electronDensities.createAttribute("zMin", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    Attribute zMax(electronDensities.createAttribute("zMax", PredType::NATIVE_DOUBLE, H5S_SCALAR));

    nZPoints.write(PredType::NATIVE_INT, &nGridPoints);
    zMin.write(PredType::NATIVE_DOUBLE, &minGridPoint);
    zMax.write(PredType::NATIVE_DOUBLE, &maxGridPoint);

    int k = 0;
    for(int j = 0; j < signed(density.n_cols); j++){
        for(int i = 0; i < signed(density.n_rows); i++){
            hsize_t dim[3] = {density(i,j).n_cols, density(i,j).n_rows, density(i,j).n_slices};
            DataSpace space(3, dim);
            stringstream orbital;
            orbital << "orbital" << setw(4) << setfill('0')  << k;
            DataSet dataset(electronDensities.createDataSet(orbital.str(), PredType::NATIVE_DOUBLE, space));
            dataset.write(density(i,j).memptr(), PredType::NATIVE_DOUBLE);
            k++;
        }
    }
}


void OutputManager::saveDipoleMoment(const double& dipoleMoment)
{

    Attribute energyAttribute(m_dataset->createAttribute("dipole Moment", PredType::NATIVE_DOUBLE, H5S_SCALAR));
    energyAttribute.write(PredType::NATIVE_DOUBLE, &dipoleMoment);

}

void OutputManager::closeOutput()
{
    m_output->flush(H5F_SCOPE_GLOBAL);
    m_output->close();
}
