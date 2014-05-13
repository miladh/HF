#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <libconfig.h++>
#include <hf.h>

using namespace arma;
using namespace std;
using namespace libconfig;
using namespace hf;

enum solverMethod {
    rhf, uhf
};

ElectronicSystem *setupSystem(string name);

int main(int argc, char **argv)
{

    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    boost::mpi::timer timer;
    timer.restart();

    //read config file---------------------------------------------------------------
    Config cfg;
    for(int p = 0; p < world.size(); p++) {

        world.barrier();
        if(p != world.rank()) {
            continue;
        }
        cfg.readFile("../../../hf/apps/default/defaultConfig.cfg");
    }
    const Setting & root = cfg.getRoot();


    //Setup system--------------------------------------------------------------------
    string chemicalSystem = root["chemicalSystem"]["name"];
    const Setting &atomsMeta = root["chemicalSystem"]["atoms"];
    vector<Atom *> atoms;

    for(int i = 0; i < atomsMeta.getLength(); i++){
        const Setting &atomMeta = atomsMeta[i];

        string basisFile;
        stringstream basisFilePath;
        rowvec position = zeros<rowvec>(3);

        atomMeta.lookupValue("basis",basisFile);
        basisFilePath << "infiles/turbomole/"<< basisFile;

        const Setting &pos =  atomMeta["position"];
        for(int i =0; i < 3; i++){
            position[i] = pos[i];
        }

        atoms.push_back(new Atom(basisFilePath.str(), position));
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);

//    ElectronicSystem *system = setupSystem("H2");

    //setup solver--------------------------------------------------------------------
    int solverMethod = root["solverSettings"]["method"];
    HFsolver* solver;
    string method;

    switch (solverMethod) {
    case rhf:
        method = "rhf";
        solver = new RHF(system);
        break;

    case uhf:
        method = "uhf";
        solver = new UHF(system);
        break;
    }



    int maxNumOfIteration = root["solverSettings"]["maxNumOfIteration"];
    double dampingFactor = root["solverSettings"]["dampingFactor"];

    solver->setDampingFactor(dampingFactor);
    solver->setMaxNumOfIteration(maxNumOfIteration);

    int useDIISprocdure = root["solverSettings"]["DIISprocedureSettings"]["useDIISprocedure"];
    if(useDIISprocdure){
        int nTerms = root["solverSettings"]["DIISprocedureSettings"]["nTerms"];
        int iterationLimit = root["solverSettings"]["DIISprocedureSettings"]["iterationLimit"];
        solver->useDIISprocedure(nTerms,iterationLimit);
    }



    //run solver--------------------------------------------------------------------
    if(world.rank()==0){
        cout << "---------------------------Hartree-Fock------------------------------"  << endl;
        cout << "system:    " << chemicalSystem << endl;
        cout << "method:    " << method << endl;
    }

    solver->runSolver();
    double laps = timer.elapsed();


    //Analyzer------------------------------------------------------------------------------
    Analyser analyser(&cfg, system,solver);
    analyser.runAnalysis();



    //Save config file-------------------------------------------------------------------------
    if(world.rank() == 0 && int(root["analysisSettings"]["saveResults"])){
        string outputFilePath = root["analysisSettings"]["outputFilePath"];
        stringstream copyCommand;
        copyCommand << "cp ../../../hf/apps/default/defaultConfig.cfg" << " " << outputFilePath;
        const char* command = new char[sizeof(copyCommand)];
        command = (copyCommand.str()).c_str();

        int status = std::system( command );
        if(!status){
            cout << "Config file successfully copied!" << endl;
        }else{
            cerr << "Config file cannot be copied!" << endl;
        }
    }


    if(world.rank()==0){
        cout << setprecision(3)
             << "Total elapsed time:  " <<  timer.elapsed()  << "s" << endl
             << " - Computation time: " << laps << "s" << endl
             << " - Analysis time:    " << timer.elapsed() - laps << "s" << endl;
    }



    return 0;

}




ElectronicSystem* setupSystem(string name)
{
    vector<Atom *> atoms;
    vector<rowvec3> atomsPos;

    if(name =="H2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -2.69, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { 2.69, 0, 0 }));

    }else if(name =="Li2"){
        atomsPos.push_back({-2.5255, 0.0, 0.0});
        atomsPos.push_back({ 2.5255, 0.0, 0.0});
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", {-2.5255, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_3_basis_3-21G.tm", { 2.5255, 0.0, 0.0}));

    }else if(name =="O2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", {-1.14, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", { 1.14, 0.0, 0.0}));

    }else if(name =="H2O"){
        double D = 1.797;
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_STO-3G.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {D, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", { -D*cos((180-104.45) *M_PI/180.0),
                                                                              D*sin((180-104.45) *M_PI/180.0), 0.0}));
    }else if(name =="CO2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", {-2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_3-21G.tm", { 2.185, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0, 0.0, 0.0}));

    }else if(name =="CO"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", {2.132, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_6-31Gd.tm", {0.0, 0.0, 0.0}));

    }else if(name =="CH4"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_6-31Gd.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {2.043/sqrt(3), 2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {-2.043/sqrt(3), -2.043/sqrt(3), 2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {2.043/sqrt(3), -2.043/sqrt(3), -2.043/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {-2.043/sqrt(3), 2.043/sqrt(3), -2.043/sqrt(3)}));

    }else if(name =="SiO4"){
        double D = 4.9;
        double T = 2.0*D /sqrt(3);
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {D/sqrt(3), D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), -D/sqrt(3), D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{D/sqrt(3), -D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-D/sqrt(3), D/sqrt(3), -D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_14_basis_3-21G.tm", {T, -T, -T}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T - D/sqrt(3), -T -D/sqrt(3), -T - D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm",{T + D/sqrt(3), -T -D/sqrt(3), -T + D/sqrt(3)}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {T + D/sqrt(3), -T +D/sqrt(3), -T -D/sqrt(3)}));

    }else if(name =="Fe2S2"){
        double FeFe = 4.556129358;
        double SS   = 6.908838214;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0,  FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31G.tm", {0.0, -FeFe*0.5, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {-SS*0.5, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_6-31G.tm", {SS*0.5, 0.0, 0.0}));

    }else if(name =="FeS2"){
        double FeS = 3.802128689;
        double SFeS = 114.0;
        atoms.push_back(new Atom("infiles/turbomole/atom_26_basis_6-31Gd.tm", {0.0, 0.0 , 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm", {FeS, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_16_basis_3-21Gd.tm",  {-FeS*cos((180.0-SFeS) *M_PI/180.0),
                                                                                FeS*sin((180.0-SFeS) *M_PI/180.0), 0.0}));

    }else if(name =="B2H6"){
        atoms.push_back(new Atom("infiles/turbomole/atom_35_basis_STO-3G.tm", {4.15523001913855,	-2.62734596571781,	-2.1668423433}));
        atoms.push_back(new Atom("infiles/turbomole/atom_35_basis_STO-3G.tm", {-1.58249411094126,	-2.71417845917264,	-1.5931496711}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {1.32268313023405,	-0.488422131801799,	-0.463250906}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {1.26964420979805,	-4.596418355359,	-3.131788174}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-4.18664215946605,	-4.97783286330769,	-2.7807851776}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-4.55222698082542,	-1.23494785798249,	-0.3039594963}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {7.50426593053485,	-0.843550494551715,	-1.3734189282}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {7.27201346643206,	-5.02418270176317,	-4.0695587434}));

    }else if(name =="H3O"){
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", {0. ,0., 0.12144702}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {0.,1.759450169,-0.32385935}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {1.523729528,-0.879726029,-0.32385935}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31Gds.tm", {-1.523729528,-0.879726029,-0.32385935}));

    }else if(name =="CH3"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.0        , 0.0,  0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.000000,    2.037899403 ,    0.000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {1.764871793 ,   -1.018949702 ,    0.000000}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { -1.764871793,    -1.018949702,     0.000000}));
    }

    else if(name =="benzene"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        ,  4.68463073,  0. }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,   0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023,  0.}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023,  0.}));


//        double R = 6.666953288;
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.         ,  4.68463073, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,  R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023, R}));
//        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023, R}));

    }else if(name =="benzeneT"){
        double R = 0.0;
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,0.000000,1.059035 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,-1.206008,1.757674 +R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,-1.207177,3.151591 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,0.000000,3.848575 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,1.207177,3.151591 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.000000,1.206008,1.757674 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,0.000000,-0.021580 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,-2.141639,1.214422 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,-2.143566,3.692995 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,0.000000,4.930150 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,2.143566,3.692995 + R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {0.000000,2.141639,1.214422 + R}));

        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {-1.394063,0.000000,-2.454152}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", { -0.697047,1.207238,-2.454628}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.697047,1.207238,-2.454628}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {1.394063,0.000000,-2.454152}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {0.697047,-1.207238,-2.454628}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_STO-3G.tm", {-0.697047,-1.207238,-2.454628}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-2.475399,0.000000,-2.450322}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-1.238232,2.143565,-2.453676}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {1.238232,2.143565,-2.453676}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {2.475399,0.000000,-2.450322}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {1.238232,-2.143565,-2.453676}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_STO-3G.tm", {-1.238232,-2.143565,-2.453676}));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);
    return system;

}

