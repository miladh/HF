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
    clock_t begin = clock();

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

//    ElectronicSystem *system = new ElectronicSystem();
//    system->addAtoms(atoms);

    ElectronicSystem *system = setupSystem("H2");

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


    //Analyzer--------------------------------------------------------------------
    Analyser analyser(system,solver);
    int saveEnergies = root["analysisSettings"]["saveEnergies"];
    int atomicPartialCharge = root["analysisSettings"]["atomicPartialCharge"];
    int chargeDensity = root["analysisSettings"]["chargeDensity"];
    int calculateElectrostaticPotential = root["analysisSettings"]["calculateElectrostaticPotential"];

    if(saveEnergies){
        analyser.saveEnergies();
    }
    if(atomicPartialCharge){
        analyser.atomicPartialCharge();
    }
    if(chargeDensity){
        analyser.calculateChargeDensity();
    }
    if(calculateElectrostaticPotential){
        analyser.calculateElectrostaticPotential();
    }
    analyser.dipoleMoment();



    clock_t end = clock();
    if(world.rank()==0){
        cout << "Total elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;
    }

    return 0;

}




ElectronicSystem* setupSystem(string name)
{
    vector<Atom *> atoms;
    vector<rowvec3> atomsPos;

    if(name =="H2"){
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { -0.69, 0, 0 }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_4-31G.tm", { 0.69, 0, 0 }));

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
        atoms.push_back(new Atom("infiles/turbomole/atom_8_basis_6-31Gds.tm", { 0.0, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", {D, 0.0, 0.0}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_6-31G.tm", { -D*cos((180-104.45) *M_PI/180.0),
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

    }else if(name =="benzene"){
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        ,  4.68463073,  0.         }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,   0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023,  0.        }));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023,  0.        }));


        double R = 6.666953288;
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        ,  2.63805748, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872,  1.31902874, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 2.28467872, -1.31902874, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", { 0.        , -2.63805748, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872, -1.31902874, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_6_basis_3-21G.tm", {-2.28467872,  1.31902874, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.         ,  4.68463073, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 ,  2.34326023, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", { 4.0572417 , -2.34326023, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {0.        , -4.68463073,  R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 , -2.34326023, R}));
        atoms.push_back(new Atom("infiles/turbomole/atom_1_basis_3-21G.tm", {-4.0572417 ,  2.34326023, R}));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }

    ElectronicSystem *system = new ElectronicSystem();
    system->addAtoms(atoms);
    return system;

}

