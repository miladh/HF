#include <iostream>
#include <armadillo>
#include <mpi.h>

#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;

System* setupSystem(string name);
void angstromToau(vector<rowvec3> &corePos);
void sampleConfigurations(System* system, HFsolver* solver);
void convertToCartesian(mat &pos, const double& r, const double& t);
int main(int argc, char **argv)
{

    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    int nProcs = MPI::COMM_WORLD.Get_size();


    clock_t begin = clock();

    /********************************************************************************/

    //options:
    string method = "rhf";
    string chemicalSystem = "CO2";

    if(rank==0){
        cout << "---------------------Configuration Sampler-----------------------"  << endl;
        cout << "system:    " << chemicalSystem << endl;
        cout << "method:    " << method << endl;
    }


    //Setup system:
    System *system = setupSystem(chemicalSystem);


    //Choose method:
    HFsolver* solver;
    if(method == "rhf"){
        solver = new RHF(system, rank, nProcs);
    }else if(method == "uhf"){
        solver = new UHF(system, rank, nProcs);
    }else{
        cerr << "unknown method!" << endl;
        exit(0);
    }


    sampleConfigurations(system, solver);

    /********************************************************************************/
    clock_t end = clock();
    if(rank==0){
        cout << "Total elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;
    }

    MPI::Finalize();
    return 0;

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

    pos.row(1) = r1;
    pos.row(2) = r2;
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

        core.push_back(new BasisSet("infiles/turbomole/O_6-31G_d"));
        core.push_back(new BasisSet("infiles/turbomole/H_6-31G_ds"));
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

void sampleConfigurations(System* system, HFsolver* solver)
{
    int nCores = system->getNumOfCores();
    mat pos = zeros(nCores, 3);
    vec bondLength = linspace(1.8,2.2,1e2);
    vec bondAngle = linspace(acos(-1)/4.0, acos(-1),1e2);
    mat data = zeros(bondLength.n_elem * bondAngle.n_elem, 3);

    int i = 0;
    for(double r: bondLength){
        for(double t: bondAngle){
            if(MPI::COMM_WORLD.Get_rank()==0){
                cout << "Bond length:    " << r << endl;
                cout << "Bond angle:    "  << t << endl;
            }
            convertToCartesian(pos, r, t);
            for(int core = 0; core < nCores; core++){
                system->m_basisSet.at(core)->setCorePosition(pos.row(core));
            }
            solver->runSolver();
            data(i,0) = solver->getEnergy();
            data(i,1) =  r;
            data(i,2) = t;
            i++;
        }

    }

    if(MPI::COMM_WORLD.Get_rank()==0){
           data.save("/home/milad/kurs/qmd/param/param.bin",raw_binary);
    }
}






























