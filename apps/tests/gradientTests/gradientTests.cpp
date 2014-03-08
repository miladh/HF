#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
using namespace hf;

//TEST(energyGradient_H2)
//{

//    //Initializing the system
//    int nElectrons = 2;
//    BasisSet *basisCoreA = new BasisSet("infiles/turbomole/H_3-21G");
//    BasisSet *basisCoreB = new BasisSet("infiles/turbomole/H_3-21G");

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

//    BOMD BOSolver(system,0,1);


//    //Domain
//    vec bondLength = linspace(0.2, 0.9, 10);
//    vec gradient   = 0*bondLength;
//    vec numericalGradient = 0*bondLength;



//    //Calculate gradient analytically
//    for(uint x = 0; x < bondLength.n_elem; x++){
//        rowvec X = {bondLength(x) , 0 ,0 };

//        system->m_basisSet.at(0)->setCorePosition(X * -0.5);
//        system->m_basisSet.at(1)->setCorePosition(X *  0.5);
//        BOSolver.solveSingleStep();

//        gradient(x) = BOSolver.getEnergyGradient()(0);
//    }


//    //Calculate gradient numerically
//    for(uint x = 0; x < bondLength.n_elem; x++){
//        rowvec X = {bondLength(x) , 0 ,0 };
//        rowvec dx = {h , 0 ,0 };

//        system->m_basisSet.at(0)->setCorePosition((X-dx) * -0.5);
//        system->m_basisSet.at(1)->setCorePosition((X-dx) *  0.5);
//        BOSolver.solveSingleStep();
//        double Ep = BOSolver.getEnergy();

//        system->m_basisSet.at(0)->setCorePosition((X+dx) * -0.5);
//        system->m_basisSet.at(1)->setCorePosition((X+dx) *  0.5);
//        BOSolver.solveSingleStep();
//        double En = BOSolver.getEnergy();

//        numericalGradient(x) = (En - Ep) / (2 * h);
//    }


//    for(uint i = 0; i < numericalGradient.n_elem; i++){
//        CHECK_CLOSE(numericalGradient(i), gradient(i), 1e-6);
////        cout << setprecision(14) << "[" << bondLength(i) << "," << numericalGradient(i) <<"," << gradient(i) << "]," <<endl;

//    }

//}


//TEST(energy_H2)
//{
//    //Initializing the system
//    int nElectrons = 2;
//    BasisSet *basisCoreA = new BasisSet("infiles/turbomole/H_3-21G");
//    BasisSet *basisCoreB = new BasisSet("infiles/turbomole/H_3-21G");

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

//    BOMD BOSolver(system,0,1);


//    //Domain
//    vec bondLength = linspace(0.2, 10.0, 1000);

//    //Calculate gradient analytically
//    for(uint x = 0; x < bondLength.n_elem; x++){
//        rowvec X = {bondLength(x) , 0 ,0 };

//        system->m_basisSet.at(0)->setCorePosition(X * -0.5);
//        system->m_basisSet.at(1)->setCorePosition(X *  0.5);
//        BOSolver.solveSingleStep();

//        cout << setprecision(14) << "[" << bondLength(x) << "," << BOSolver.getEnergy() << "]," <<endl;
//    }

//}



