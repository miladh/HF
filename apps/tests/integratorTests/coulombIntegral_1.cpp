#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;


SUITE(DEVELOPMENT){
TEST(GTOnuclearAttractionIntegral)
{
    /*
     * test case:   electron-nucleus integral
     * max angular momentum:    2
     *
     * source:
     *      numerical integration
     *      python scripts
     * */

    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec posC = {2.3,0.9,3.2};


    PrimitiveGTO primitiveA(0.2, 1.0);
    PrimitiveGTO primitiveB(0.3, 1.0);
    PrimitiveGTO primitiveC(0.0, 0.0);
    primitiveA.setCenter(posA);
    primitiveB.setCenter(posB);
    integrator.setCorePositionC(posC);


    integrator.setMaxAngularMomentum(2);


    bool oneParticle = true;
    bool twoParticle = false;
    bool kinetic = false;

    //Build E-cube for Lmax = 2: max power in x,y and z is 2!
    primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    integrator.updateHermiteCoefficients(oneParticle, twoParticle, kinetic);



    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.788948987251e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.971203468743e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.024071525839e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.727033700014e-03, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.134361291529e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.921666495443e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.185957751329e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.105746642202e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(9.596523510045e-03, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.388444040338e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-9.204700718547e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.019226499202e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-5.276399683274e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.927318225377e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-6.299787002318e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-9.718105595370e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.037280861539e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.312309453843e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.202910466576e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,1}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.073449397904e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.319499900436e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.435114042344e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.536931448007e+00, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.067865861209e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.033153544029e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.524622701603e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.703919381580e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.292169308884e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.162963233448e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,0,2}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.390872806284e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.637350724302e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.139721853567e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.213713540367e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.136233458775e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.306634838230e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.697963144320e-04, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.907709578263e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.932098923684e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.414126847830e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.842342257899e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.356912442721e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.187325132374e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.128036711345e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.083519267815e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.544897767601e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.894393296797e-04, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.132616876012e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.386353605834e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.911527548440e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,1,1}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.229240947242e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.609849112824e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(9.032384498820e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.625292648498e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.939589748931e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.913397190183e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.878262296370e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.131065513841e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.052929737663e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.267402937768e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({0,2,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.289710831960e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.786414780578e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-9.322262110550e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.671155215998e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.222106053447e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.972830178046e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.026352276293e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.576450345257e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-3.945885129414e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-4.918734877201e-03, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.459437143524e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.263894353489e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.735756798558e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.071773603054e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.115384967585e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.802219023191e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.350111738398e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.197526800952e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.145639876363e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.638641395940e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,0,1}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.278875254192e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.185667243163e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.417205698627e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.560020091608e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-2.926456829930e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-7.176178735649e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-5.223967979758e-04, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(9.269318129877e-03, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.337071697343e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(-1.203714316117e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({1,1,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.401501778682e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(7.889586550718e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.935977010010e-01, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,2});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(5.534914541236e-01, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(2.563391673303e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,1,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.217850538435e-02, integrator.nuclearAttractionIntegral(), 1e-5);


    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,2,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(8.419480232293e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(1.481688684288e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,0,1});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(3.878852644576e-02, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,1,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(4.176920693786e-03, integrator.nuclearAttractionIntegral(), 1e-5);

    primitiveA.setPowers({2,0,0}); primitiveB.setPowers({2,0,0});
    integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
    CHECK_CLOSE(6.422210627967e-02, integrator.nuclearAttractionIntegral(), 1e-5);

}
}
