#include <unittest++/UnitTest++.h>
#include <armadillo>
#include <iostream>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<bomd/bomd.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/splitValence/h_321g.h>

using namespace std;
using namespace arma;
void writeToFile(const vec X, const vec E, const vec F);

//TEST(energyGradient)
//{
//    vec bondLength= linspace(0.5, 3.0, 1e3);
//    vec gradient  = 0*bondLength;
//    vec energy    = 0*bondLength;

//    int nElectrons = 2;

//    BasisSet *basisCoreA = new H_321G;
//    BasisSet *basisCoreB = new H_321G;

//    rowvec coreCharges = {1 , 1};
//    rowvec coreMass = {1 , 1};

//    int maxAngularMomentum = basisCoreA->getAngularMomentum()+1;
//    basisCoreA->setCorePosition(zeros<rowvec>(3));
//    basisCoreA->setCoreCharge(coreCharges(0));
//    basisCoreA->setCoreMass(coreMass(0));

//    basisCoreB->setCorePosition(zeros<rowvec>(3));
//    basisCoreB->setCoreCharge(coreCharges(1));
//    basisCoreB->setCoreMass(coreMass(1));

//    System *system = new System(nElectrons, maxAngularMomentum);
//    system->addBasisSet(basisCoreA);
//    system->addBasisSet(basisCoreB);

//    BOMD BOSolver(system);

//    for(uint x = 0; x < bondLength.n_elem; x++){
//        rowvec X = {bondLength(x) , 0 ,0 };
//        system->m_basisSet.at(0)->setCorePosition(X * -0.5);
//        system->m_basisSet.at(1)->setCorePosition(X *  0.5);


//        BOSolver.solveSingleStep();
//        energy(x) = BOSolver.getEnergy();
//        gradient(x) = BOSolver.getEnergyGradient()(0);
//        cout << "[" << bondLength(x) << "," << energy(x) <<"," << gradient(x) << "]," <<endl;
//    }

//    writeToFile(bondLength, energy, gradient);

//}

void writeToFile(const vec X, const vec E, const vec F){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/energyForce.dat";
    myfile.open(outName.str().c_str(),ios::binary);

    for(uint i=0;  i < X.n_elem; i++){
        myfile  << X(i) << "  " << E(i) << "  " << F(i) << endl;
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}
