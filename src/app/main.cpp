#include <iostream>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/h_321g.h>
#include<basisSet/li_321g.h>
#include<basisSet/o_321g.h>


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


    int m_case = 5;

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
        A = {-0.5, 0.0, 0.0};
        B = {0.5, 0.0, 0.0};
        coreCharges = {1 , 1};
        basisCoreA  = new H_321G;
        basisCoreB  = new H_321G;

    }else if(m_case==3){
        //Lithium molecule
        nElectrons = 6;
        A = {-2.5, 0.0, 0.0};
        B = { 2.5, 0.0, 0.0};
        coreCharges = {3 , 3};
        basisCoreA = new Li_321G;
        basisCoreB = new Li_321G;

    }else if(m_case==4){
        //Oxygen molecule
        nElectrons = 16;
        A = {-1.14, 0.0, 0.0};
        B = { 1.14, 0.0, 0.0};
        coreCharges = {8 , 8};
        basisCoreA = new O_321G;
        basisCoreB = new O_321G;

    }else if(m_case==5){
        //Water molecule
        nElectrons = 10;
        A = {-2.0, 4.0, 0.0};
        B = { 2.0, 4.0, 0.0};
        C = { 0.0, 0.0, 0.0};
        coreCharges = {1 , 1, 8};
        basisCoreA = new H_321G;
        basisCoreB = new H_321G;
        basisCoreC = new O_321G;

    }

    /********************************************************/

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCorePosition(A);
    basisCoreB->setCorePosition(B);

    if(m_case == 5){
        maxAngularMomentum = basisCoreC->getAngularMomentum();
        basisCoreC->setCorePosition(C);
    }

    System system(nElectrons, maxAngularMomentum, coreCharges);
    system.addBasisSet(basisCoreA);
    system.addBasisSet(basisCoreB);

    if(m_case == 5){
        system.addBasisSet(basisCoreC);
    }

    HFsolver solver(system);
    solver.runSolver();


    double end_time = time(NULL);
    cout << "-------------------------------"  << endl;
    cout << "Elapsed time: " << end_time-start_time << "s" << endl;


    return 0;
}


