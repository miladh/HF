#include <iostream>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<cpmd/cpmd.h>
#include<bomd/bomd.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/splitValence/h_321g.h>
#include<basisSet/splitValence/h_431g.h>
#include<basisSet/splitValence/li_321g.h>
#include<basisSet/splitValence/o_321g.h>
#include<basisSet/splitValence/o_431g.h>
#include<basisSet/splitValence/c_321g.h>


using namespace arma;
using namespace std;


int main()
{
    clock_t begin = clock();
    int nElectrons;
    rowvec coreCharges,coreMass, A, B, C;

    BasisSet *basisCoreA;
    BasisSet *basisCoreB;
    BasisSet *basisCoreC;


    int m_case = 2;
    int dynamic = 0;
    int cpmd = 0;

    if(m_case == 1){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.5, 0.0, 0.0};
        B = {0.5, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA  = new H_QuadZeta;
        basisCoreB  = new H_QuadZeta;

    }else if(m_case==2){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-2.25, 0.0, 0.0};
        B = {2.25, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA  = new H_321G;
        basisCoreB  = new H_321G;

    }else if(m_case==3){
        //Hydrogen molecule
        nElectrons = 2;
        A = {-0.69, 0.0, 0.0};
        B = {0.69, 0.0, 0.0};
        coreCharges = {1 , 1};
        coreMass = {1 , 1};
        basisCoreA  = new H_431G;
        basisCoreB  = new H_431G;

    }else if(m_case==4){
        //Lithium molecule
        nElectrons = 6;
        A = {-2.5255, 0.0, 0.0};
        B = { 2.5255, 0.0, 0.0};
        coreCharges = {3 , 3};
        coreMass = {7 , 7};
        basisCoreA = new Li_321G;
        basisCoreB = new Li_321G;

    }else if(m_case==5){
        //Oxygen molecule
        nElectrons = 16;
        A = {-1.14, 0.0, 0.0};
        B = { 1.14, 0.0, 0.0};
        coreMass = {16 , 16};
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
        coreMass = {1 , 1, 16};
        basisCoreA = new H_431G;
        basisCoreB = new H_431G;
        basisCoreC = new O_431G;


    }else if(m_case==7){
        //Carbon dioxide
        nElectrons = 22;
        A = {-2.2, 0.0, 0.0};
        B = { 2.2, 0.0, 0.0};
        C = { 0.0, 0.0, 0.0};
        coreCharges = {8 , 8, 6};
        coreMass = {16 , 16, 12};
        basisCoreA = new O_321G;
        basisCoreB = new O_321G;
        basisCoreC = new C_321G;

    }


    /********************************************************/

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCoreCharge(coreCharges(0));
    basisCoreA->setCoreMass(coreMass(0));
    basisCoreA->setCorePosition(A);

    basisCoreB->setCorePosition(B);
    basisCoreB->setCoreCharge(coreCharges(1));
    basisCoreB->setCoreMass(coreMass(1));


    if(m_case == 6 ||m_case == 7 ){
        maxAngularMomentum = basisCoreC->getAngularMomentum();
        basisCoreC->setCorePosition(C);
        basisCoreC->setCoreCharge(coreCharges(2));
        basisCoreC->setCoreMass(coreMass(2));

    }


    if(dynamic){
        maxAngularMomentum +=1;
    }

    System *system = new System(nElectrons, maxAngularMomentum);
    system->addBasisSet(basisCoreA);
    system->addBasisSet(basisCoreB);

    if(m_case == 6 || m_case == 7 ){
        system->addBasisSet(basisCoreC);
    }

    if(dynamic){
        if(cpmd){
          CPMD cpSolver(system);
          cpSolver.runDynamics();

        }else{
            BOMD boSolver(system);
            boSolver.runDynamics();
        }
    }else{
        HFsolver solver(system);
        solver.runSolver();
    }


    sleep(1.2);
    clock_t end = clock();
    cout << "-------------------------------"  << endl;
    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << "s" << endl;


    return 0;
}


