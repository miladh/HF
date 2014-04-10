#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;


SUITE(DEVELOPMENT){
TEST(electronRepulsionIntegral)
{
    /*
     * test case:   electron-electron integral
     * max angular momentum:    2
     *
     * source:
     *      numerical integration
     *      python scripts
     * */

    Integrator integrator(2);

    const rowvec posA = {1.2,2.3,3.4};
    const rowvec posB = {-1.3,1.4,-2.4};
    const rowvec posC = {2.3,0.9,3.2};
    const rowvec posD = {5.0,1.9,1.2};

    PrimitiveGTO primitiveA(0.2, 1.0);
    PrimitiveGTO primitiveB(0.3, 1.0);
    PrimitiveGTO primitiveC(0.4, 1.0);
    PrimitiveGTO primitiveD(0.1, 1.0);
    primitiveA.setCenter(&posA);
    primitiveB.setCenter(&posB);
    primitiveC.setCenter(&posC);
    primitiveD.setCenter(&posD);


    //Build E-cube for Lmax = 2: max power in x,y and z is 2!
    primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
    primitiveC.setPowers({2,2,2}); primitiveD.setPowers({2,2,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.setPrimitiveC(primitiveC);integrator.setPrimitiveD(primitiveD);
    integrator.updateOverlapHermiteCoefficients();
    integrator.updateElectronRepulsionHermiteCoefficients();


    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,0});
    primitiveC.setPowers({0,0,0}); primitiveD.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.setPrimitiveC(primitiveC);integrator.setPrimitiveD(primitiveD);
    CHECK_CLOSE(0.1624848293, integrator.electronRepulsionIntegral(),1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,0});
    primitiveC.setPowers({0,0,0}); primitiveD.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.setPrimitiveC(primitiveC);integrator.setPrimitiveD(primitiveD);
    CHECK_CLOSE(0.2667434786, integrator.electronRepulsionIntegral(),1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,0});
    primitiveC.setPowers({0,2,0}); primitiveD.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.setPrimitiveC(primitiveC);integrator.setPrimitiveD(primitiveD);
    CHECK_CLOSE(0.2681206721, integrator.electronRepulsionIntegral(),1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,0,0});
    primitiveC.setPowers({0,2,0}); primitiveD.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.setPrimitiveC(primitiveC);integrator.setPrimitiveD(primitiveD);
    CHECK_CLOSE(-0.8837373783, integrator.electronRepulsionIntegral(),1e-5);
   }
}
