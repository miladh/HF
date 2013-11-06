#include <iostream>
#include <armadillo>

#include<system/system.h>
#include<hfSolver/hfsolver.h>
#include<basisSet/h_quadzeta.h>
#include<basisSet/h_321g.h>
#include<basisSet/o_321g.h>
#include<basisSet/h_sto6.h>


using namespace arma;
using namespace std;


int main()
{
    double start_time = time(NULL);

//    H_QuadZeta *basisCoreA  = new H_QuadZeta;
//    H_QuadZeta *basisCoreB  = new H_QuadZeta;

//    H_321G *basisCoreA = new H_321G;
//    H_321G *basisCoreB = new H_321G;

    O_321G *basisCoreA = new O_321G;
    O_321G *basisCoreB = new O_321G;

//    H_STO6 *basisCoreA = new H_STO6;
//    H_STO6 *basisCoreB = new H_STO6;

    int nElectrons = 16;
    rowvec coreCharges = {8 , 8};
    rowvec A = {-1.14 , 0, 0};
    rowvec B = { 1.14 , 0, 0};

    /********************************************************/

    int maxAngularMomentum = basisCoreA->getAngularMomentum();
    basisCoreA->setCorePosition(A);
    basisCoreB->setCorePosition(B);

    System system(nElectrons, maxAngularMomentum, coreCharges);
    system.addBasisSet(basisCoreA);
    system.addBasisSet(basisCoreB);

    HFsolver solver(system);
    solver.runSolver();


    double end_time = time(NULL);
    cout << "-------------------------------"  << endl;
    cout << "Elapsed time: " << end_time-start_time << "s" << endl;


    return 0;
}


