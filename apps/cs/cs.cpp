#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <H5Cpp.h>
#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;

int main(int argc, char* argv[])
{

    boost::mpi::environment env;
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();


    // Read in output file, abort if there are too few command-line arguments
    if( argc < 2 ){
        cerr << "Bad Usage: input file is not specified!" << endl;
        exit(1);
    }
    string inFileName = argv[1];

    // MPI----------------------------------------------------------------------
    vector<int> stateIDs;
    if(world.rank() == 0) {
        H5::H5File inFile(inFileName, H5F_ACC_RDONLY );
        H5::Group statesGroup(inFile.openGroup("/states"));
        int nStates = statesGroup.getNumObjs();
        vector<int> allStates;
        allStates.reserve(nStates);
        for(int i = 0; i < nStates; i++) {
            allStates.push_back(i);
        }
        random_shuffle(allStates.begin(), allStates.end());
        for(int p = 0; p < world.size(); p++) {
            vector<int> states;
            for(int i = BLOCK_LOW(p, signed(world.size()), signed(allStates.size()));
                i <= BLOCK_HIGH(p, signed(world.size()), signed(allStates.size()));
                i++) {

                states.push_back(allStates.at(i));

            }
            if(p == 0) {
                stateIDs = states;
            } else {
                world.send(p, 0, states);
            }
        }
        inFile.close();
    }
    else {
        world.recv(0, 0, stateIDs);
    }

    cout << "Rank " << world.rank() << ": I have " << stateIDs.size() << " states to compute" << endl;
    world.barrier();
    //------------------------------------------------------------------------

    //Create output file
    stringstream outFileName;
    outFileName << "/home/milad/kurs/qmd/tmpData/uhf_" << world.rank() << ".h5";
    H5::H5File outputFile(outFileName.str(), H5F_ACC_TRUNC);


    //Read input file and copy metadata, states and attributes
    for(int p = 0; p < world.size(); p++) {

        world.barrier();
        if(p != world.rank()) {
            continue;
        }

        //Read input file:
        H5::H5File inputFile(inFileName, H5F_ACC_RDONLY);


        //Loop over input file objects:
        for(int objectIndex = 0; objectIndex < int(inputFile.getNumObjs()); objectIndex++) {
            string objectName = inputFile.getObjnameByIdx(objectIndex);

            if(objectName == "states") {
                H5::Group statesGroupIn(inputFile.openGroup(objectName));
                H5::Group statesGroupOut(outputFile.createGroup(objectName));

                for(int stateID : stateIDs){
                    string stateName = statesGroupIn.getObjnameByIdx(stateID);
                    H5Ocopy(statesGroupIn.getId(), stateName.c_str(), statesGroupOut.getId(), stateName.c_str(),
                            H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
                }

            }else{
                H5Ocopy(inputFile.getId(), objectName.c_str(), outputFile.getId(), objectName.c_str(),
                        H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
            }
        }
        inputFile.close();
    }



    //-------------------------------------------------------------------------------------------------------

    //Create additional dataset elements

    struct AtomData {
        double x;
        double y;
        double z;
        double q;
    };

    H5::CompType atomCompound( sizeof(AtomData) );
    atomCompound.insertMember("x", HOFFSET(AtomData, x), H5::PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("y", HOFFSET(AtomData, y), H5::PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("z", HOFFSET(AtomData, z), H5::PredType::NATIVE_DOUBLE);
    atomCompound.insertMember("q", HOFFSET(AtomData, q), H5::PredType::NATIVE_DOUBLE);


    struct AtomMetaData {
        int    type;
        char basisName[64];
    };
    H5::StrType string_type(H5::PredType::C_S1, 64);
    H5::CompType atomMetaCompound( sizeof(AtomMetaData) );
    atomMetaCompound.insertMember( "type", HOFFSET(AtomMetaData, type), H5::PredType::NATIVE_INT);
    atomMetaCompound.insertMember( "basisName", HOFFSET(AtomMetaData, basisName), string_type);


    //-------------------------------------------------------------------------------------------------------

    //Read meta data
    H5::Group rootGroup(outputFile.openGroup("/"));
    H5::DataSet atomMeta(rootGroup.openDataSet("atomMeta"));
    hsize_t dims[1];
    atomMeta.getSpace().getSimpleExtentDims(dims);
    int nAtoms = dims[0];

    AtomMetaData atomMetaData[nAtoms];
    atomMeta.read(atomMetaData, atomMetaCompound);


    //Read state data
    H5::Group statesGroup(outputFile.openGroup("/states"));

    for(int stateID = 0; stateID < signed(statesGroup.getNumObjs()); stateID++) {
        string stateName = statesGroup.getObjnameByIdx(stateID);
        H5::DataSet stateDataSet(statesGroup.openDataSet(stateName));
        hsize_t dims2[1];
        stateDataSet.getSpace().getSimpleExtentDims(dims2);
        int nAtoms2 = dims2[0];

        if(nAtoms != nAtoms2) {
            throw std::logic_error("Number of atoms in state data not equal to the one given in meta data!");
        }

        AtomData *atoms = new AtomData[nAtoms2];
        stateDataSet.read(atoms, atomCompound);

        ElectronicSystem system;
        vector<Atom *> atomList;
        for(int i = 0; i < nAtoms2; i++) {
            stringstream basisFile;
            basisFile << "infiles/turbomole/atom_" << atomMetaData[i].type
                      << "_basis_" << atomMetaData[i].basisName << ".tm";
            atomList.push_back(new Atom(basisFile.str(), { atoms[i].x, atoms[i].y, atoms[i].z}));
        }
        system.addAtoms(atomList);


        //Choose method:
        HFsolver* solver;
        string method = "uhf";

        if(method == "rhf"){
            solver = new RHF(&system);
        }else if(method == "uhf"){
            solver = new UHF(&system);
        }else{
            cerr << "unknown method!" << endl;
            exit(0);
        }
        solver->runSolver();

        //-------------------------------------------------------------------------------------------------------
        Analyser analyser(&system, solver);
        analyser.computeAtomicPartialCharge();

        for(int i = 0; i < nAtoms2; i++) {
            Atom* atom = atomList.at(i);
            atoms[i].q = atom->corePartialCharge();
            stateDataSet.write(atoms, atomCompound);
        }

        double energy = solver->energy();
        H5::Attribute energyAttribute(stateDataSet.createAttribute("energy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR));
        energyAttribute.write(H5::PredType::NATIVE_DOUBLE, &energy);





        outputFile.flush(H5F_SCOPE_GLOBAL);

        delete atoms;
    }


    world.barrier();
    if(world.rank() == 0) {
        cout << "Total time "  << fixed << setprecision(2) << timer.elapsed() << " s" << endl;
    }
    outputFile.close();
    return 0;

}























