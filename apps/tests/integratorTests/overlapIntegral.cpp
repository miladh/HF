#include <unittest++/UnitTest++.h>
#include <hf.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;
using namespace hf;
SUITE(DEVELOPMENT){
    TEST(GTOoverlapIntegral)
    {
        /*
     * test case:   overlap integral S
     * max angular momentum:    2
     *
     * source:
     *      numerical integration
     *      python scripts
     * */


        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};

        PrimitiveGTO primitiveA(0.2, 1.0);
        PrimitiveGTO primitiveB(0.3, 1.0);

        primitiveA.setCenter(&posA);
        primitiveB.setCenter(&posB);

        Integrator integrator(2);

        //Build E-cube for Lmax = 2: max power in x,y and z is 2!
        primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        integrator.updateOverlapHermiteCoefficients();


        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.191723635809e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.764798835076e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(7.606056933184e-01, integrator.overlapIntegral(), 1e-5);


        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.290205088911e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(9.953275806273e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.346171019009e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.191723635809e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.764798835076e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.290205088911e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.383447271617e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-8.429776310255e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-2.093948045733e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-3.034719471692e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.684675146152e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-8.429776310255e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(), 1e-5);


        primitiveA.setPowers({0,0,1}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-8.294396505227e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.562397355490e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.795322214215e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(6.361589630418e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(5.624630479765e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.006315997117e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.764884052762e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.562397355490e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.795322214215e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(5.624630479765e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,0,2}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.124794710981e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.107270743920e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(9.600525610073e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.227321941537e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.311086675171e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(9.600525610073e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.287061526673e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.130731944696e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-6.791027795542e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.562581629595e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,1,1}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.478974112823e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.539230248010e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.571014175384e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(9.823983134901e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-7.329386373895e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.700417638744e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.195477460565e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.539230248010e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.571014175384e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-7.329386373895e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({0,2,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.078460496021e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.787585453713e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.140908539978e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(), 1e-5);


        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-2.019256528514e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-5.958618179043e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.382399417538e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-2.145102544455e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.191723635809e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(6.220797378920e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.264466446538e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.140922068599e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(7.027012719229e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.073599126307e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.214888155128e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(7.464956854705e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,0,1}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.147198252614e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(9.652961450049e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(6.160906115879e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.440078841511e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-1.966630012757e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.217653816683e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(7.464956854705e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(-4.800262805037e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({1,1,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(6.435307633366e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.873101816378e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(8.985596213996e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.471968503285e+00, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.394316653896e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(3.234814637039e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(4.375055811780e-01, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.979309089521e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(6.911997087690e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(1.072551272228e-02, integrator.overlapIntegral(), 1e-5);

        primitiveA.setPowers({2,0,0}); primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_CLOSE(2.979309089521e-01, integrator.overlapIntegral(), 1e-5);
    }




    /*************************************************************************************
 *
 * Derivative of overlap integral
 *
 *************************************************************************************/


    TEST(GTOoverlapIntegral_derivative)
    {

        Integrator integrator(2);

        rowvec posA = {1.2,2.3,3.4};
        rowvec posB = {-1.3,1.4,-2.4};
        rowvec Sab  = {0, 0, 0};


        PrimitiveGTO primitiveA(0.2, 1.0);
        PrimitiveGTO primitiveB(0.3, 1.0);
        primitiveA.setCenter(&posA);
        primitiveB.setCenter(&posB);

        //Build E-cube for Lmax = 2: max power in x,y and z is 2!
        primitiveA.setPowers({2,2,2}); primitiveB.setPowers({2,2,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);

        integrator.updateOverlapHermiteCoefficients();
        integrator.updateOverlapHermiteCoefficientsGD();



        Sab  = {-7.150341814851e-02,-2.574123053346e-02,-1.658879301045e-01};
        primitiveA.setPowers({0,0,0}); primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(Sab, integrator.QDerivativeOverlapIntegral(),3, 1e-5);


        Sab  = {-7.150341814851e-02,-2.574123053346e-02,-1.658879301045e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.658879301045e-01,-5.971965483764e-02,-3.371910524102e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.563634159911e-01,-1.642908297568e-01,-8.375792182932e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.574123053346e-02,3.840210244029e-02,-5.971965483764e-02};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-5.971965483764e-02,8.909287766148e-02,-1.213887788677e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-8.077026114056e-02,5.244346700684e-03,-1.873870058461e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.383447271617e-02,-2.574123053346e-02,-1.658879301045e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-5.529597670152e-02,-5.971965483764e-02,-3.371910524102e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-8.580410177821e-03,3.840210244029e-02,-5.971965483764e-02};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.766894543234e-02,-5.148246106693e-02,-3.317758602091e-01};
        primitiveA.setPowers({0,0,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.488318951568e-01,8.957948225645e-02,5.057865786153e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {5.057865786153e-01,1.820831683015e-01,8.416490021784e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.256368827440e+00,4.522927778783e-01,1.784030158849e+00};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {8.957948225646e-02,-1.336393164922e-01,1.820831683015e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.820831683015e-01,-2.716411118217e-01,3.029936407842e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.810805087691e-01,-1.825032651838e-02,5.713365192039e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {8.294396505227e-02,8.957948225646e-02,5.057865786153e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.685955262051e-01,1.820831683015e-01,8.416490021784e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.658879301045e-01,1.791589645129e-01,1.011573157231e+00};
        primitiveA.setPowers({0,0,1});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-9.374384132942e-01,-3.374778287859e-01,-1.677193328529e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.677193328529e+00,-6.037895982704e-01,-2.254556422760e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-3.816953778251e+00,-1.374103360170e+00,-4.106337339291e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-3.374778287859e-01,5.034669238332e-01,-6.037895982704e-01};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.037895982704e-01,9.007646303086e-01,-8.116403121938e-01};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.058930431657e+00,6.875548298465e-02,-1.894557583906e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-3.124794710981e-01,-3.374778287859e-01,-1.677193328529e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-5.590644428430e-01,-6.037895982704e-01,-2.254556422760e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.124926095953e-01,5.034669238332e-01,-6.037895982704e-01};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.249589421962e-01,-6.749556575719e-01,-3.354386657058e+00};
        primitiveA.setPowers({0,0,2});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {3.861184580020e-02,-5.760315366044e-02,8.957948225646e-02};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {8.957948225646e-02,-1.336393164922e-01,1.820831683015e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.464362446352e-01,-3.676463679224e-01,4.522927778783e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-5.760315366044e-02,-7.221959638469e-02,-1.336393164922e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.336393164922e-01,-1.675494636125e-01,-2.716411118217e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-7.866520051026e-03,-6.798003478342e-03,-1.825032651838e-02};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.287061526673e-02,-5.760315366044e-02,8.957948225646e-02};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.920105122015e-02,-7.221959638469e-02,-1.336393164922e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.574123053346e-02,-1.152063073209e-01,1.791589645129e-01};
        primitiveA.setPowers({0,1,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.343692233847e-01,2.004589747383e-01,-2.731247524523e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.731247524523e-01,4.074616677325e-01,-4.544904611764e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.784391668175e-01,1.012130727385e+00,-9.633762857784e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.004589747383e-01,2.513241954187e-01,4.074616677325e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {4.074616677325e-01,5.108525369867e-01,6.780324361549e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.737548977757e-02,2.365705210462e-02,5.564461623294e-02};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-9.104158415075e-02,4.074616677325e-01,-4.544904611764e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {6.681965824611e-02,2.513241954187e-01,4.074616677325e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-8.957948225645e-02,4.009179494767e-01,-5.462495049045e-01};
        primitiveA.setPowers({0,1,1});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-9.235381488062e-02,4.397631824337e-02,-2.142608505230e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.142608505230e-01,1.020250583246e-01,-4.355159632930e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-5.894389880941e-01,2.806744535565e-01,-1.081817318347e+00};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {4.397631824337e-02,-3.780562283286e-02,1.020250583246e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.020250583246e-01,-8.770904497223e-02,2.073805898175e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.917286476339e-01,-1.433904442414e-01,-4.448104625106e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-3.078460496021e-02,4.397631824337e-02,-2.142608505230e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-7.142028350768e-02,1.020250583246e-01,-4.355159632930e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.465877274779e-02,-3.780562283286e-02,1.020250583246e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.156920992041e-02,8.795263648674e-02,-4.285217010461e-01};
        primitiveA.setPowers({0,2,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {3.575170907426e-02,3.861184580020e-02,2.488318951568e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {8.294396505227e-02,8.957948225646e-02,5.057865786153e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.281817079955e-01,2.464362446352e-01,1.256368827440e+00};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.287061526673e-02,-5.760315366044e-02,8.957948225646e-02};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {4.038513057028e-02,-7.866520051026e-03,2.810805087691e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.072551272228e-01,1.287061526673e-02,8.294396505227e-02};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.488318951568e-01,2.985982741882e-02,1.685955262051e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-3.861184580020e-02,-1.920105122015e-02,2.985982741882e-02};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.191723635809e-01,2.574123053346e-02,1.658879301045e-01};
        primitiveA.setPowers({1,0,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.244159475784e-01,-1.343692233847e-01,-7.586798679230e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.528932893077e-01,-2.731247524523e-01,-1.262473503268e+00};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.281844137199e-01,-6.784391668175e-01,-2.676045238273e+00};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-9.104158415076e-02,4.074616677325e-01,-4.544904611764e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.405402543846e-01,2.737548977757e-02,-8.570047788058e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {3.732478427352e-01,-4.478974112823e-02,-2.528932893077e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {7.586798679230e-01,-9.104158415075e-02,-4.208245010892e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.343692233847e-01,6.681965824611e-02,-9.104158415076e-02};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {4.147198252614e-01,-8.957948225645e-02,-5.057865786153e-01};
        primitiveA.setPowers({1,0,1});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.930592290010e-02,8.640473049066e-02,-1.343692233847e-01};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.232181223176e-01,5.514695518836e-01,-6.784391668175e-01};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.880157683022e-02,1.083293945770e-01,2.004589747383e-01};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {6.681965824611e-02,2.513241954187e-01,4.074616677325e-01};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {3.933260025514e-03,1.019700521751e-02,2.737548977757e-02};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {5.791776870029e-02,2.880157683022e-02,-4.478974112823e-02};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {1.343692233847e-01,6.681965824611e-02,-9.104158415076e-02};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-8.640473049066e-02,3.610979819234e-02,6.681965824611e-02};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {6.435307633366e-02,5.760315366044e-02,-8.957948225645e-02};
        primitiveA.setPowers({1,1,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.787585453713e-02,-8.365899923376e-02,-5.391357728398e-01};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-4.147198252614e-02,-1.940888782223e-01,-1.095870920333e+00};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.140908539978e-01,-5.339451967095e-01,-2.722132459453e+00};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,0,2});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-6.435307633366e-03,1.248068329310e-01,-1.940888782223e-01};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.492991370941e-02,2.895518523998e-01,-3.945135313199e-01};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,1,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-2.019256528514e-02,1.704412677723e-02,-6.090077689998e-01};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({0,2,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {2.085516362665e-01,-6.435307633366e-03,-4.147198252614e-02};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {4.838397961383e-01,-1.492991370941e-02,-8.429776310255e-02};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,0,1});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {7.507858905594e-02,9.600525610074e-03,-1.492991370941e-02};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({1,1,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

        Sab  = {-1.191723635809e-02,-6.435307633366e-02,-4.147198252614e-01};
        primitiveA.setPowers({2,0,0});primitiveB.setPowers({2,0,0});
        integrator.setPrimitiveA(primitiveA);integrator.setPrimitiveB(primitiveB);
        CHECK_ARRAY_CLOSE(integrator.QDerivativeOverlapIntegral(),Sab,3, 1e-5);

    }
}
