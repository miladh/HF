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


System* setupSystem(string name);
void angstromToau(vector<rowvec3> &corePos);
void sampleConfigurations(System* system, HFsolver* solver);
void convertToCartesian(mat &pos, const double& r, const double& t);


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
            for(int i = blockLow(p, world.size(), allStates.size());
                i <= blockHigh(p, world.size(), allStates.size());
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
    outFileName << "output_" << world.rank() << ".h5";
    H5::H5File outputFile(outFileName.str(), H5F_ACC_TRUNC);


    //Read input file and copy metadata, states and attributes
    for(int p = 0; p < world.size(); p++) {

        world.barrier();
        if(p != world.rank()) {
            continue;
        }

        //Read input file:
        H5::H5File inputFile(argv[1], H5F_ACC_RDONLY);


        //Loop over input file objects:
        for(int objectIndex = 0; objectIndex < int(inputFile.getNumObjs()); objectIndex++) {
            string objectName = inputFile.getObjnameByIdx(objectIndex);

            if(objectName == "states") {
                H5::Group statesGroup(inputFile.openGroup(objectName));
                H5::Group statesGroupOut(outputFile.createGroup(objectName));

                for(int stateID : stateIDs){
                    string stateName = statesGroup.getObjnameByIdx(stateID);
                    H5Ocopy(statesGroup.getId(), stateName.c_str(), statesGroupOut.getId(), stateName.c_str(),
                            H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
                }

            }else{
                H5Ocopy(inputFile.getId(), objectName.c_str(), outputFile.getId(), objectName.c_str(),
                        H5Pcreate(H5P_OBJECT_COPY), H5P_DEFAULT);
            }
        }
        inputFile.close();
    }





    /********************************************************************************/

    //options:
    //    string method = "rhf";
    //    string chemicalSystem = "CO2";


    //    //Setup system:
    //    System *system = setupSystem(chemicalSystem);

    //    if(rank==0){
    //        cout << "---------------------Configuration Sampler-----------------------"  << endl;
    //        cout << "system:    " << chemicalSystem << endl;
    //        cout << "method:    " << method << endl;
    //    }

    //    //Choose method:
    //    HFsolver* solver;
    //    if(method == "rhf"){
    //        solver = new RHF(system, rank, nProcs);
    //    }else if(method == "uhf"){
    //        solver = new UHF(system, rank, nProcs);
    //    }else{
    //        cerr << "unknown method!" << endl;
    //        exit(0);
    //    }


    //    sampleConfigurations(system, solver);

    //    /********************************************************************************/
    if(world.rank()==0){
        cout << "Total elapsed time: "<< timer.elapsed() << "s" << endl;
    }


    outputFile.close();

    return 0;

}

void sampleConfigurations(System* system, HFsolver* solver)
{
    int nCores = system->getNumOfCores();
    int Nr = 5e1;
    double rMin = 2.0;
    double rMax = 2.4;

    int Nt = 5e1;
    double tMin = 3.0*acos(-1)/4.0;
    double tMax = acos(-1);

    vec bondLength = linspace(rMin, rMax, Nr);
    vec bondAngle = linspace(tMin, tMax, Nt);
    mat pos = zeros(nCores, 3);
    mat data = zeros(bondLength.n_elem * bondAngle.n_elem, 3);

    int i = 0;
    for(double r: bondLength){
        for(double t: bondAngle){
            data(i,1) =  r;
            data(i,2) =  t;
            i++;
        }

    }


    int myRank = 0;
    int nProcs = 1;
    // MPI----------------------------------------------------------------------
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    vector<int> myGridIndices;

    int node = 0;
    int s = 0;
    for (int p = 0; p < Nr; p++) {
        if (myRank == node){
            myGridIndices.push_back(p);
        }
        s++;
        if(s >= blockSize(node, nProcs, Nr)){
            s = 0;
            node++;
        }
    }

    cout << "Rank: " << myRank << " - Number of grid points: " << myGridIndices.size() << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    //---------------------------------------------------------------------------

    i = 0;
    for(int r: myGridIndices){
        i = Nt * r;
        for(double t: bondAngle){
            cout << "Bond length:    " <<  bondLength[r] << endl;
            cout << "Bond angle:    "  << t << endl;

            convertToCartesian(pos, bondLength[r], t);
            for(int core = 0; core < nCores; core++){
                system->m_basisSet.at(core)->setCorePosition(pos.row(core));
            }
            solver->runSolver();
            data(i, 0) = solver->getEnergy();
            i++;
        }

    }

    stringstream fileName;
    fileName << "/home/milad/kurs/qmd/param/id" << myRank << "_param.bin";
    data.save(fileName.str(),raw_binary);
}










System* setupSystem(string name)
{
    int nElectrons;
    rowvec coreCharges,coreMass;
    vector<BasisSet*> core;
    vector<rowvec3> corePos;

    if(name =="CO2"){
        nElectrons = 22;
        coreCharges = {8 , 8, 6};
        coreMass = {16 , 16, 12};
        corePos.push_back({-2.2, 0.0, 0.0});
        corePos.push_back({ 2.2, 0.0, 0.0});
        corePos.push_back({ 0.0, 0.0, 0.0});
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/C_3-21G"));

    }else if(name =="H2O"){
        nElectrons = 10;
        coreCharges = {8 , 1, 1};
        coreMass = {16 , 1, 1};
        nElectrons = 10;
        corePos.push_back({ 0.0, 0.0, 0.0});
        corePos.push_back({1.797, 0.0, 0.0});
        corePos.push_back({ -1.797*cos((180-104.45) *M_PI/180.0),
                            1.797*sin((180-104.45) *M_PI/180.0), 0.0});

        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/H_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/H_6-31G_ds"));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    int maxAngularMomentum = core[0]->getAngularMomentum();

    System *system = new System(nElectrons, maxAngularMomentum);

    for (uint i = 0; i < core.size(); i++){
        core[i]->setCorePosition(corePos[i]);
        core[i]->setCoreCharge(coreCharges(i));
        core[i]->setCoreMass(coreMass(i));
        system->addBasisSet(core[i]);
    }

    return system;
}



void angstromToau(vector<rowvec3>& corePos)
{
    for(rowvec3& pos: corePos){
        pos *= 1.889725989;
    }
}

void convertToCartesian(mat &pos, const double& r, const double& t)
{
    rowvec r1 = {-r * sin(t*0.5), r * cos(t*0.5), 0};
    rowvec r2 = { r * sin(t*0.5), r * cos(t*0.5), 0};

    pos.row(0) = r1;
    pos.row(1) = r2;
}


























