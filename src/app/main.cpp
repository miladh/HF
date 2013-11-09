#include <iostream>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/splitValence/h_321g.h>
#include<basisSet/splitValence/h_431g.h>
#include<basisSet/splitValence/li_321g.h>
#include<basisSet/splitValence/o_321g.h>
#include<basisSet/splitValence/o_431g.h>


using namespace arma;
using namespace std;


int main()
{
    double start_time = time(NULL);
    int nElectrons;
    rowvec coreCharges, A, B, C;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;
    BasisSet *basisCoreC;


    int m_case = 6;

    if(m_case == 1){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.5, 0.0, 0.0};
        B = {0.5, 0.0, 0.0};
        coreCharges = {1 , 1};
        basisCoreA  = new H_QuadZeta;
        basisCoreB  = new H_QuadZeta;

    }else if(m_case==2){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.7, 0.0, 0.0};
        B = {0.7, 0.0, 0.0};
        coreCharges = {1 , 1};
        basisCoreA  = new H_321G;
        basisCoreB  = new H_321G;

    }else if(m_case==3){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.69, 0.0, 0.0};
        B = {0.69, 0.0, 0.0};
        coreCharges = {1 , 1};
        basisCoreA  = new H_431G;
        basisCoreB  = new H_431G;

    }else if(m_case==4){
        //Lithium molecule
        nElectrons = 6;
        A = {-2.5255, 0.0, 0.0};
        B = { 2.5255, 0.0, 0.0};
        coreCharges = {3 , 3};
        basisCoreA = new Li_321G;
        basisCoreB = new Li_321G;

    }else if(m_case==5){
        //Oxygen molecule
        nElectrons = 16;
        A = {-1.14, 0.0, 0.0};
        B = { 1.14, 0.0, 0.0};
        coreCharges = {8 , 8};
        basisCoreA = new O_431G;
        basisCoreB = new O_431G;

    }else if(m_case==6){
        double x = 1.797*cos((180-104.45) *M_PI/180.0);
        double y = 1.797*sin((180-104.45) *M_PI/180.0);
        //Water molecule
        nElectrons = 10;
        A = {1.797, 0.0, 0.0};
        B = { -x, y, 0.0};
        C = { 0.0, 0.0, 0.0};
        coreCharges = {1 , 1, 8};
        basisCoreA = new H_431G;
        basisCoreB = new H_431G;
        basisCoreC = new O_431G;

    }

    /********************************************************/

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCorePosition(A);
    basisCoreB->setCorePosition(B);

    if(m_case == 6){
        maxAngularMomentum = basisCoreC->getAngularMomentum();
        basisCoreC->setCorePosition(C);
    }

    System system(nElectrons, maxAngularMomentum, coreCharges);
    system.addBasisSet(basisCoreA);
    system.addBasisSet(basisCoreB);

    if(m_case == 6){
        system.addBasisSet(basisCoreC);
    }

    HFsolver solver(system);
    solver.runSolver();


    double end_time = time(NULL);
    cout << "-------------------------------"  << endl;
    cout << "Elapsed time: " << end_time-start_time << "s" << endl;


    return 0;
}



//vec X = linspace(0.5, 4, 20);
//for(uint i = 0; i < X.n_elem; i++)
//{

//    double x = X[i]*cos((180-104.45) *M_PI/180.0);
//    double y = X[i]*sin((180-104.45) *M_PI/180.0);
//    //Water molecule
//    nElectrons = 10;
//    A = {X[i], 0.0, 0.0};
//    B = { -x, y, 0.0};
//    C = { 0.0, 0.0, 0.0};
//    coreCharges = {1 , 1, 8};


//    basisCoreA = new H_431G;
//    basisCoreB = new H_431G;
//    basisCoreC = new O_431G;

//    basisCoreA->setCorePosition(A);
//    basisCoreB->setCorePosition(B);
//    basisCoreC->setCorePosition(C);
//    int maxAngularMomentum = basisCoreC->getAngularMomentum();
//    System system(nElectrons, maxAngularMomentum, coreCharges);
//    system.addBasisSet(basisCoreA);
//    system.addBasisSet(basisCoreB);
//    system.addBasisSet(basisCoreC);
//    HFsolver solver(system);
//    solver.runSolver();

//}

