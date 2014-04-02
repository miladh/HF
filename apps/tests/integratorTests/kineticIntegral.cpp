#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;
SUITE(DEVELOPMENT){
TEST(GTOkineticIntegral)
{
    /*
     * test case:   kinetic integral K
     * max angular momentum:    2
     *
     * source:
     *      numerical integration
     *      python scripts
     * */

    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};

    PrimitiveGTO primitiveA(0.2, 1.0);
    PrimitiveGTO primitiveB(0.3, 1.0);
    primitiveA.setCenter(posA);
    primitiveB.setCenter(posB);

    integrator.setMaxAngularMomentum(2);


    bool oneParticle = true;
    bool twoParticle = false;
    bool kinetic = true;

    //Build E-cube for Lmax = 2: max power in x,y and z is 2!
    primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.updateHermiteCoefficients(oneParticle, twoParticle, kinetic);


    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-9.678702680582e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.581907301477e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.289130997844e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.454683743671e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.306080091812e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.209847292592e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-6.818565954642e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-9.183555810588e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.425034522333e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.554388972658e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.372860952215e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.514020826620e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.256885959785e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.959120137717e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.767148243772e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.085946543142e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.377533371588e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.908745121589e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.418618463595e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-6.190706910076e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.907944812827e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.109128468812e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-8.787431724837e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.368298260464e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-8.520922147503e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.440953256899e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.800828501289e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.711681426600e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.381742282583e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.682025615507e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.959120137717e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.903857602951e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-6.049047312582e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-8.688217105502e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.938468694693e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.137551783499e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.744921166165e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.304752788336e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.438680206576e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.650722365658e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(9.550382489858e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.303232565825e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.887481505265e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.840890648636e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.274267732432e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.013966668721e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.711001832636e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.512318573258e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.651535687867e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-5.921817515938e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.606516375488e-03, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.708290533872e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.565122085492e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.142903313736e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.794492285775e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.598401092187e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.532083467153e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.022784893196e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.377533371588e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.195516000820e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.137551783499e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.330149372044e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-8.809221115897e-03, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-5.361497900979e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-8.319565708416e-03, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.818565954642e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.066300057382e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.363117682384e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.652884024961e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.274267732432e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.942436762098e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.042246851468e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.789282215612e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.686798511658e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.377533371588e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.206327675248e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.930388273683e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.617381749247e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.013966668720e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.379694245101e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.247934856262e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.686798511658e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.583728636544e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.137551783498e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.287537353407e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.150543568546e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-5.213314751101e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.888774502915e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.578401717871e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.962760567236e-01, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.704641488661e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.295888952647e-02, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.562586305833e-03, integrator.kineticIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.468392469657e-01, integrator.kineticIntegral(), 1e-5);
}



/*************************************************************************************
 *
 * Derivative of kinetic integral
 *
 *************************************************************************************/



TEST(GTOkineticIntegral_derivative)
{

    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec Tab  = {0, 0, 0};
    PrimitiveGTO primitiveA(0.2, 1.0);
    PrimitiveGTO primitiveB(0.3, 1.0);
    primitiveA.setCenter(posA);
    primitiveB.setCenter(posB);

    integrator.setMaxAngularMomentum(2);


    bool oneParticle = true;
    bool twoParticle = false;
    bool kinetic = true;

    //Build E-cube for Lmax = 2: max power in x,y and z is 2!
    primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.updateHermiteCoefficients(oneParticle, twoParticle, kinetic);
    integrator.updateHermiteCoefficients_derivative(oneParticle, twoParticle, kinetic);


    Tab  = {4.091139572785e-02,1.472810246203e-02,9.491443808862e-02};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.510133486353e-02,1.983648055087e-02,1.005608330648e-01};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {8.782064003279e-02,3.161543041181e-02,1.302754383914e-01};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {8.550207133995e-03,-2.419618925033e-02,1.983648055087e-02};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.503763389837e-03,-3.475286842201e-02,7.068592975089e-03};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.320597488176e-02,7.753874778772e-03,1.234378617257e-01};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-3.523688446359e-03,8.550207133996e-03,5.510133486353e-02};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.144599160392e-02,5.503763389837e-03,1.963498048636e-02};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-3.327826283366e-03,-1.497968466466e-02,5.503763389837e-03};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.727426381857e-02,2.121901115334e-02,1.367447385438e-01};
    primitiveA.setPowers({0,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.265200229529e-02,-2.975472082630e-02,-1.508412495972e-01};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.945247072954e-02,-1.060288946263e-02,4.187293763460e-02};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.061153609984e-01,3.820152995943e-02,3.732782385369e-01};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.255645084755e-03,5.212930263301e-02,-1.060288946264e-02};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.309707092973e-02,3.154992602106e-02,8.779273133667e-02};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.176974704839e-01,-3.136356259454e-02,-2.198521566410e-01};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.216898740587e-02,-8.255645084755e-03,-2.945247072955e-02};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.115712886245e-01,3.309707092973e-02,2.438686981574e-01};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092972e-02};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-5.510133486352e-02,-3.084400733054e-02,-1.398307940375e-01};
    primitiveA.setPowers({0,0,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.464571954140e-01,5.272459034902e-02,1.744766887699e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.280497100773e-01,-8.209789562783e-02,-6.856183990998e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-9.826166149089e-01,-3.537419813672e-01,-1.763373539942e+00};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.827008855960e-02,-1.078153621574e-01,-8.209789562774e-02};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.270073992127e-01,7.031047633496e-02,-4.416162986024e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.571120252525e-01,1.302537805037e-01,3.611049576180e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.761661540526e-01,-2.827008855960e-02,-2.280497100771e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-4.785429688727e-01,-2.270073992127e-01,-1.226711940562e+00};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-9.041804176183e-02,1.301669956255e-02,-2.270073992127e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {9.763813027596e-02,-2.543724513444e-03,-1.877484875895e-01};
    primitiveA.setPowers({0,0,2});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.282531070099e-02,3.629428387549e-02,-2.975472082630e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.255645084755e-03,5.212930263301e-02,-1.060288946263e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.172155309473e-02,9.204039914690e-02,3.820152995943e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.246952699699e-02,2.518944398691e-02,5.212930263301e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.005586667488e-02,1.822763878263e-02,3.154992602106e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.351877698041e-02,1.837984583950e-02,-3.136356259454e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {4.991739425050e-03,2.246952699699e-02,-8.255645084754e-03};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092973e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-6.334914546175e-03,7.856740854583e-03,2.005586667488e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.550207133995e-03,5.415555857965e-02,-3.084400733054e-02};
    primitiveA.setPowers({0,1,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.238346762713e-02,-7.819395394952e-02,1.590433419396e-02};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-4.964560639459e-02,-4.732488903159e-02,-1.316890970050e-01};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.201276949753e-01,2.686939080730e-02,-4.327805573967e-01};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-3.008380001232e-02,-2.734145817395e-02,-4.732488903160e-02};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.046591122420e-02,8.436186990016e-02,1.571926473900e-01};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.361546143843e-02,-5.828417101632e-02,1.108690935568e-01};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.209847605341e-02,5.046591122421e-02,-2.407668076873e-01};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.808222059976e-02,3.297634872656e-02,5.046591122419e-02};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {8.255645084751e-03,-9.224103598278e-02,-5.559125239680e-02};
    primitiveA.setPowers({0,1,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {6.857419882415e-02,-9.639098252967e-04,1.590921412720e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.076695371465e-01,2.224974320323e-02,2.040771936304e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.138436938137e-01,1.004268359624e-01,3.525579867858e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {9.590406553117e-03,2.675467803669e-03,2.224974320323e-02};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {4.673575720113e-02,-1.484308548881e-02,1.020333523379e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.078924496974e-01,-9.834808072654e-04,2.503104832979e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {6.931507033685e-04,9.590406553116e-03,1.076695371465e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.553275841004e-02,4.673575720113e-02,9.955336244010e-02};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.375111856278e-02,-6.397881676221e-03,4.673575720115e-02};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {4.571613254943e-02,1.214460218731e-02,2.496208103767e-01};
    primitiveA.setPowers({0,2,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {5.285532669538e-03,-1.282531070099e-02,-8.265200229529e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.216898740587e-02,-8.255645084755e-03,-2.945247072954e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.203805097404e-01,1.172155309473e-02,1.061153609984e-01};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {4.991739425050e-03,2.246952699699e-02,-8.255645084755e-03};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092973e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.474306569721e-03,-1.351877698041e-02,-1.176974704839e-01};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {6.136709359178e-02,4.991739425050e-03,3.216898740587e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {8.265200229529e-02,1.874719404663e-02,1.115712886245e-01};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.282531070099e-02,-6.334914546176e-03,1.874719404663e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.670319847949e-02,-8.550207133995e-03,-5.510133486352e-02};
    primitiveA.setPowers({1,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-4.825348110881e-02,1.238346762713e-02,4.417870609432e-02};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.673569329367e-01,-4.964560639459e-02,-3.658030472361e-01};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-5.053504583775e-01,-2.201276949754e-01,-1.202168214991e+00};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-8.209847605341e-02,5.046591122420e-02,-2.407668076873e-01};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-4.234024791492e-02,5.361546143843e-02,1.240970880481e-01};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.239780034429e-01,-2.812079106995e-02,-1.673569329367e-01};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-4.417870609430e-02,-8.209847605341e-02,-4.249279898629e-01};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.238346762713e-02,3.808222059976e-02,-8.209847605341e-02};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {4.140562735410e-02,8.255645084759e-03,2.945247072955e-02};
    primitiveA.setPowers({1,0,1});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-7.487609137574e-03,-3.370429049548e-02,1.238346762713e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-9.457782461603e-02,-5.707906268285e-03,-2.201276949753e-01};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {9.502371819266e-03,-1.178511128187e-02,-3.008380001232e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.808222059976e-02,3.297634872655e-02,5.046591122421e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {9.591335708574e-03,-2.512248750705e-02,5.361546143843e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.923796605149e-02,9.502371819265e-03,-2.812079106995e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-1.238346762712e-02,3.808222059975e-02,-8.209847605341e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {3.370429049548e-02,2.207068427119e-02,3.808222059976e-02};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {6.425011141151e-03,-2.246952699699e-02,8.255645084759e-03};
    primitiveA.setPowers({1,1,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.022784893196e-02,2.933264701749e-02,1.890326141127e-01};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.377533371588e-02,2.147041030723e-02,8.404514972528e-02};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.195516000820e-02,-1.553924858651e-02,-1.796618623990e-01};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {2.137551783499e-03,-5.098810130656e-02,2.147041030723e-02};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.375940847460e-03,-4.879995045527e-02,-6.442699361569e-02};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.330149372044e-02,2.897602265551e-02,2.662544063993e-01};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {-2.065018716129e-02,2.137551783499e-03,1.377533371588e-02};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {6.821311685898e-02,1.375940847467e-03,4.908745121587e-03};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.058479399536e-02,-3.744921166164e-03,1.375940847460e-03};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

    Tab  = {1.097834880885e-01,5.947253902451e-02,3.832674737135e-01};
    primitiveA.setPowers({2,0,0});primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_ARRAY_CLOSE(integrator.kineticIntegral_derivative(),Tab,3, 1e-5);

}
}
