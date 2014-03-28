#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <H5Cpp.h>
#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;


int blockLow(int id, int np, int n) {
    return (id * n) / np;
}

int blockHigh(int id, int np, int n) {
    return blockLow(id + 1, np, n) - 1;
}

int blockSize(int id, int p, int n) {
    return blockLow(id + 1,p,n) - blockLow(id,p,n);
}


ElectronicSystem* setupSystem(string name);

int main(int argc, char* argv[])
{

//    boost::mpi::environment env;
//    boost::mpi::communicator world;
//    boost::mpi::timer timer;
//    timer.restart();


//    // Read in output file, abort if there are too few command-line arguments
//    if( argc < 2 ){
//        cerr << "Bad Usage: input file is not specified!" << endl;
//        exit(1);
//    }
//    string inFileName = argv[1];

//    // MPI----------------------------------------------------------------------

//    rowvec A = {1,1,1};
//    Atom atom("infiles/turbomole/C_3-21G",A);

//    exit(0);

//    // MPI----------------------------------------------------------------------
//    vector<int> stateIDs;
//    if(world.rank() == 0) {
//        H5::H5File inFile(inFileName, H5F_ACC_RDONLY );
//        H5::Group statesGroup(inFile.openGroup("/states"));
//        int nStates = statesGroup.getNumObjs();
//        vector<int> allStates;
//        allStates.reserve(nStates);
//        for(int i = 0; i < nStates; i++) {
//            allStates.push_back(i);
//        }
//        random_shuffle(allStates.begin(), allStates.end());
//        for(int p = 0; p < world.size(); p++) {
//            vector<int> states;
//            for(int i = blockLow(p, world.size(), allStates.size());
//                i <= blockHigh(p, world.size(), allStates.size());
//                i++) {

//                states.push_back(allStates.at(i));

//            }
//            if(p == 0) {
//                stateIDs = states;
//            } else {
//                world.send(p, 0, states);
//            }
//        }
//        inFile.close();
//    }
//    else {
//        world.recv(0, 0, stateIDs);
//    }

//    cout << "Rank " << world.rank() << ": I have " << stateIDs.size() << " states to compute" << endl;
//    world.barrier();
//    //------------------------------------------------------------------------

//    //Create output file
//    stringstream outFileName;
//    outFileName << "output_" << world.rank() << ".h5";
//    H5::H5File outputFile(outFileName.str(), H5F_ACC_TRUNC);


//    //Read input file and copy metadata, states and attributes
//    for(int p = 0; p < world.size(); p++) {

//        world.barrier();
//        if(p != world.rank()) {
//            continue;
//        }

//        //Read input file:
//        H5::H5File inputFile(argv[1], H5F_ACC_RDONLY);


//        //Loop over input file objects:
//        for(int objectIndex = 0; objectIndex < int(inputFile.getNumObjs()); objectIndex++) {
//            string objectName = inputFile.getObjnameByIdx(objectIndex);

//            if(objectName == "states") {
//                H5::Group statesGroup(inputFile.openGroup(objectName));
//                H5::Group statesGroupOut(outputFile.createGroup(objectName));

//                for(int stateID : stateIDs){
//                    string stateName = statesGroup.getObjnameByIdx(stateID);
//                    H5Ocopy(statesGroup.getId(), stateName.c_str(), statesGroupOut.getId(), stateName.c_str(),
//                            H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
//                }

//            }else{
//                H5Ocopy(inputFile.getId(), objectName.c_str(), outputFile.getId(), objectName.c_str(),
//                        H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
//            }
//        }
//        inputFile.close();
//    }



    //-------------------------------------------------------------------------------------------------------

//    //Create additional dataset elements
//    struct AtomData {
//        double partialCharge;
//    };

//    H5::CompType atomCompound( sizeof(AtomData) );
//    atomCompound.insertMember("partialCharge", HOFFSET(AtomData, partialCharge), H5::PredType::NATIVE_DOUBLE);

//    struct AtomData {
//        double x;
//        double y;
//        double z;
//        double partialCharge;
//    };

//    H5::CompType atomCompound( sizeof(AtomData) );
//    atomCompound.insertMember("x", HOFFSET(AtomData, x), H5::PredType::NATIVE_DOUBLE);
//    atomCompound.insertMember("y", HOFFSET(AtomData, y), H5::PredType::NATIVE_DOUBLE);
//    atomCompound.insertMember("z", HOFFSET(AtomData, z), H5::PredType::NATIVE_DOUBLE);


//    struct AtomMetaData {
//        int    type;
//        char basisName[64];
//    };
//    H5::StrType string_type(H5::PredType::C_S1, 64);
//    H5::CompType atomMetaCompound( sizeof(AtomMetaData) );
//    atomMetaCompound.insertMember( "type", HOFFSET(AtomMetaData, type), H5::PredType::NATIVE_INT);
//    atomMetaCompound.insertMember( "basisName", HOFFSET(AtomMetaData, basisName), string_type);
//    H5::Group rootGroup(outputFile.openGroup("/"));
//    H5::DataSet atomMeta(rootGroup.openDataSet("atomMeta"));

//    hsize_t dims[1];
//    atomMeta.getSpace().getSimpleExtentDims(dims);
//    int nAtoms = dims[0];

//    AtomMetaData atomMetaData[nAtoms];
//    atomMeta.read(atomMetaData, atomMetaCompound);

//    H5::Group statesGroup(outputFile.openGroup("/states"));

//    int nTotal = statesGroup.getNumObjs();
//    int currentState = 0;
//    for(int stateID = 0; stateID < nTotal; stateID++) {

//        string stateName = statesGroup.getObjnameByIdx(stateID);
//        H5::DataSet stateDataSet(statesGroup.openDataSet(stateName));
//        hsize_t dims2[1];
//        stateDataSet.getSpace().getSimpleExtentDims(dims2);
//        int nAtoms2 = dims2[0];


//        AtomData *atoms = new AtomData[nAtoms2];
//        stateDataSet.read(atoms, atomCompound);


//        int nElectrons = 0;
//        int maxAngularMomentum  = 0;
//        vector<int> coreCharges;
//        vector<int> coreMass;
//        vector<BasisSet*> core;
//        vector<rowvec3> corePos;

//        for(int i = 0; i < nAtoms2; i++) {
//            stringstream basisFile;
//            nElectrons += atomMetaData[i].type;
//            basisFile << "infiles/turbomole/atom_" << atomMetaData[i].type
//                      << "_basis_" << atomMetaData[i].basisName << ".tm";
//            string fileName = basisFile.str();


//            core.push_back(new BasisSet(fileName));
//            corePos.push_back({ atoms[i].x, atoms[i].y, atoms[i].z});
//            coreCharges.push_back(nElectrons);
//            coreMass.push_back(nElectrons);
//            maxAngularMomentum = max(maxAngularMomentum, core[i]->getAngularMomentum());


//        }

//        System *system = new System(nElectrons, maxAngularMomentum);

//        for (uint i = 0; i < core.size(); i++){
//            core[i]->setCorePosition(corePos[i]);
//            core[i]->setCoreCharge(coreCharges[i]);
//            core[i]->setCoreMass(coreMass[i]);
//            system->addBasisSet(core[i]);
//        }


//        //Choose and run solver
//        HFsolver* solver;
//        string method = "rhf";
//        if(method == "rhf"){
//            solver = new RHF(system, world.rank(), world.size());
//        }else if(method == "uhf"){
//            solver = new UHF(system,  world.rank(), world.size());
//        }else{
//            cerr << "unknown method!" << endl;
//            exit(0);
//        }

//        solver->runSolver();
//        double energy = solver->getEnergy();

//        H5::Attribute energyAttribute(stateDataSet.createAttribute("energy", H5::PredType::NATIVE_DOUBLE, H5S_SCALAR));
//        energyAttribute.write(H5::PredType::NATIVE_DOUBLE, &energy);
//        outputFile.flush(H5F_SCOPE_GLOBAL);

//        currentState++;

//        delete atoms;
//    }


//    if(world.rank()==0){
//        cout << "Total elapsed time: "<< timer.elapsed() << "s" << endl;
//    }

//    outputFile.close();
//    return 0;

}

ElectronicSystem *setupSystem(string name)
{
//    int nElectrons;
//    rowvec coreCharges,coreMass;
//    vector<BasisSet*> core;
//    vector<rowvec3> corePos;

//    if(name =="CO2"){
//        nElectrons = 22;
//        coreCharges = {8 , 8, 6};
//        coreMass = {16 , 16, 12};
//        corePos.push_back({-2.2, 0.0, 0.0});
//        corePos.push_back({ 2.2, 0.0, 0.0});
//        corePos.push_back({ 0.0, 0.0, 0.0});
//        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
//        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
//        core.push_back(new BasisSet("infiles/turbomole/C_3-21G"));

//    }else if(name =="H2O"){
//        nElectrons = 10;
//        coreCharges = {8 , 1, 1};
//        coreMass = {16 , 1, 1};
//        nElectrons = 10;
//        corePos.push_back({ 0.0, 0.0, 0.0});
//        corePos.push_back({1.797, 0.0, 0.0});
//        corePos.push_back({ -1.797*cos((180-104.45) *M_PI/180.0),
//                            1.797*sin((180-104.45) *M_PI/180.0), 0.0});

//        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
//        core.push_back(new BasisSet("infiles/turbomole/H_3-21G"));
//        core.push_back(new BasisSet("infiles/turbomole/H_6-31G_ds"));

//    }else{
//        cerr << "unknown system!" << endl;
//        exit(0);
//    }

//    int maxAngularMomentum = core[0]->getAngularMomentum();

//    System *system = new System(nElectrons, maxAngularMomentum);

//    for (uint i = 0; i < core.size(); i++){
//        core[i]->setCorePosition(corePos[i]);
//        core[i]->setCoreCharge(coreCharges(i));
//        core[i]->setCoreMass(coreMass(i));
//        system->addBasisSet(core[i]);
//    }

//    return system;
}


























