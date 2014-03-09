#include <iostream>
#include <armadillo>
#include <mpi.h>

#include <hf.h>

using namespace arma;
using namespace std;
using namespace hf;

System* setupSystem(string name, int dynamic);

int main(int argc, char **argv)
{

    MPI::Init(argc, argv);
    int rank = MPI::COMM_WORLD.Get_rank();
    int nProcs = MPI::COMM_WORLD.Get_size();

    clock_t begin = clock();


    int dynamic = 0;
    System *system = setupSystem("H2",dynamic);

    if(dynamic)
    {
        BOMD boSolver(system, rank, nProcs);
        boSolver.runDynamics();

    }else{
        RHF *solver = new RHF(system, rank, nProcs);
        solver->runSolver();
    }


    clock_t end = clock();
    if(rank==0){
        cout << "-------------------------------"  << endl;
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;
    }

    MPI::Finalize();
    return 0;

}


System* setupSystem(string name, int dynamic=0)
{
    int nElectrons;
    rowvec coreCharges,coreMass;
    vector<BasisSet*> core;
    vector<rowvec3> corePos;

    if(name =="H2"){
        nElectrons = 2;
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        corePos.push_back({ -0.5, 0, 0 });
        corePos.push_back({  0.5, 0, 0 });
        core.push_back(new BasisSet("infiles/turbomole/H_Qzeta"));
        core.push_back(new BasisSet("infiles/turbomole/H_Qzeta"));

    }else if(name =="Li2"){
        nElectrons = 6;
        coreCharges = {3 , 3};
        coreMass = {7 , 7};
        corePos.push_back({-2.5255, 0.0, 0.0});
        corePos.push_back({ 2.5255, 0.0, 0.0});
        core.push_back(new BasisSet("infiles/turbomole/Li_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/Li_3-21G"));

    }else if(name =="O2"){
        nElectrons = 8;
        coreMass = {16 , 16};
        coreCharges = {8 , 8};
        corePos.push_back({-1.14, 0.0, 0.0});
        corePos.push_back({ 1.14, 0.0, 0.0});
        core.push_back(new BasisSet("infiles/turbomole/O_4-31G"));
        core.push_back(new BasisSet("infiles/turbomole/O_4-31G"));

    }else if(name =="H2O"){
        nElectrons = 10;
        coreCharges = {8 , 1, 1};
        coreMass = {16 , 1, 1};
        nElectrons = 10;
        corePos.push_back({1.797, 0.0, 0.0});
        corePos.push_back({ -1.797*cos((180-104.45) *M_PI/180.0),
                            1.797*sin((180-104.45) *M_PI/180.0), 0.0});
        corePos.push_back({ 0.0, 0.0, 0.0});

        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/H_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/H_3-21G"));

    }else if(name =="CO2"){
        nElectrons = 22;
        coreCharges = {8 , 8, 6};
        coreMass = {16 , 16, 12};
        corePos.push_back({-2.2, 0.0, 0.0});
        corePos.push_back({ 2.2, 0.0, 0.0});
        corePos.push_back({ 0.0, 0.0, 0.0});
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/C_3-21G"));


    }else if(name =="CH4"){
        nElectrons = 10;
        coreCharges = {6, 1 , 1, 1, 1};
        coreMass = {6 , 1, 1, 1, 1};
        corePos.push_back({0.0, 0.0, 0.0});
        corePos.push_back({2.043/sqrt(3), 2.043/sqrt(3), 2.043/sqrt(3)});
        corePos.push_back({-2.043/sqrt(3), -2.043/sqrt(3), 2.043/sqrt(3)});
        corePos.push_back({2.043/sqrt(3), -2.043/sqrt(3), -2.043/sqrt(3)});
        corePos.push_back({-2.043/sqrt(3), 2.043/sqrt(3), -2.043/sqrt(3)});
        core.push_back(new BasisSet("infiles/turbomole/C_4-31G"));
        core.push_back(new BasisSet("infiles/turbomole/H_4-31G"));
        core.push_back(new BasisSet("infiles/turbomole/H_4-31G"));
        core.push_back(new BasisSet("infiles/turbomole/H_4-31G"));
        core.push_back(new BasisSet("infiles/turbomole/H_4-31G"));


    }else if(name =="SiO4"){
        nElectrons = 46;
        coreCharges = {14, 8 , 8, 8, 8};
        coreMass = {28 , 16, 16, 16, 16};
        corePos.push_back({0.0, 0.0, 0.0});
        corePos.push_back({4.9/sqrt(3), 4.9/sqrt(3), 4.9/sqrt(3)});
        corePos.push_back({-4.9/sqrt(3), -4.9/sqrt(3), 4.9/sqrt(3)});
        corePos.push_back({4.9/sqrt(3), -4.9/sqrt(3), -4.9/sqrt(3)});
        corePos.push_back({-4.9/sqrt(3), 4.9/sqrt(3), -4.9/sqrt(3)});
        core.push_back(new BasisSet("infiles/turbomole/Si_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/O_3-21G"));

    }else if(name =="Fe2S2"){
        nElectrons = 84;
        coreCharges = {26, 26 , 16, 16};
        coreMass = {56 , 56, 32, 32};
        corePos.push_back({0.0, 1.0 , 0.0});
        corePos.push_back({1.0 , 0.0, 0.0});
        corePos.push_back({0.0, 0.0, 0.0});
        corePos.push_back({1.0, 1.0, 0.0});
        core.push_back(new BasisSet("infiles/turbomole/Fe_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/Fe_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/S_3-21G"));
        core.push_back(new BasisSet("infiles/turbomole/S_3-21G"));

    }else{
        cerr << "unknown system!" << endl;
        exit(0);
    }


    int maxAngularMomentum = core[0]->getAngularMomentum();
    if(dynamic){
        maxAngularMomentum += 1;
    }

    System *system = new System(nElectrons, maxAngularMomentum);

    for (uint i = 0; i < core.size(); i++){
        core[i]->setCorePosition(corePos[i]);
        core[i]->setCoreCharge(coreCharges(i));
        core[i]->setCoreMass(coreMass(i));
        system->addBasisSet(core[i]);
    }


    return system;


}

