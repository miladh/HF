#include <unittest++/UnitTest++.h>
#include <src/integrator/integrator.h>
#include <src/math/boys.h>


TEST(GTOIntegration)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);

    integrator.addPrimitives(new PrimitiveGTO(0.2, 1.0));
    integrator.addPrimitives(new PrimitiveGTO(0.3, 1.0));

    integrator.setMaxAngularMomentum(2);

    integrator.setupE();

    CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,0), 0.119172363580852, 1e-5);
    CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,1), 0.276479883507577, 1e-5);
    CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,0,2), 0.760605693318432, 1e-5);
    CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,0), 0.0429020508891068, 1e-5);
    CHECK_CLOSE(integrator.overlapIntegral(0,0,0,0,1,1), 0.0995327580627279, 1e-5);

    CHECK_CLOSE(-0.0967870268058250, integrator.kineticIntegral(0,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(-0.158190730147696, integrator.kineticIntegral(0,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(-0.0245468374367114, integrator.kineticIntegral(0,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-0.0330608009181156, integrator.kineticIntegral(0,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-0.0681856595464203, integrator.kineticIntegral(0,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-0.0918355581058773, integrator.kineticIntegral(0,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-0.0142503452233260, integrator.kineticIntegral(0,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-0.00917293898306074, integrator.kineticIntegral(0,0,0,1,1,1), 1e-5);
    CHECK_CLOSE(0.251402082662018, integrator.kineticIntegral(0,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(0.0495912013771734, integrator.kineticIntegral(0,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(0.0176714824377215, integrator.kineticIntegral(0,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(0.137753337158815, integrator.kineticIntegral(0,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(0.0490874512158928, integrator.kineticIntegral(0,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(0.0137594084745900, integrator.kineticIntegral(0,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(-0.0551617848828804, integrator.kineticIntegral(0,0,1,1,1,1), 1e-5);
    CHECK_CLOSE(-0.0604904731258242, integrator.kineticIntegral(0,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(-0.0868821710550227, integrator.kineticIntegral(0,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(0.0213755178349884, integrator.kineticIntegral(0,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(0.0137594084745909, integrator.kineticIntegral(0,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(-0.0374492116616455, integrator.kineticIntegral(0,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(-0.0334264444581330, integrator.kineticIntegral(0,1,0,1,1,1), 1e-5);
    CHECK_CLOSE(0.0788748150526478, integrator.kineticIntegral(0,1,1,0,1,1), 1e-5);
    CHECK_CLOSE(-0.0206391127118871, integrator.kineticIntegral(0,1,1,1,0,0), 1e-5);
    CHECK_CLOSE(0.0827426773243232, integrator.kineticIntegral(0,1,1,1,0,1), 1e-5);
}

TEST(BoysfactorialFunctions)
{
    Boys boys;
    CHECK_EQUAL(boys.factorial(0),1);
    CHECK_EQUAL(boys.factorial(1),1);
    CHECK_EQUAL(boys.factorial(5),120);
    CHECK_EQUAL(boys.factorial(-1),1);

    CHECK_EQUAL(boys.doubleFactorial(0),1);
    CHECK_EQUAL(boys.doubleFactorial(1),1);
    CHECK_EQUAL(boys.doubleFactorial(5),15);
    CHECK_EQUAL(boys.doubleFactorial(-1),1);
}


TEST(BoysDownwardrecursionFunction)
{
    Boys boys;
    double Fp, x;
    vec xvec, Fpvec;


    //test for Fn(0):
    x = 0;
    for(int n = 10; n > 0; n--){
        Fp = boys.downwardRecursion(x,n);
    }
    CHECK_CLOSE(Fp, 1.0, 1e-10);


    //test for F0(x) for small x:
    xvec = linspace(0.01,0.09,10);
    Fpvec = zeros(10);
    for(uint i = 0; i < 10; i++){
        for(uint n = 10; n > 0; n--){
            Fpvec[i] = boys.downwardRecursion(xvec[i],n);
        }
    }
    CHECK_CLOSE(Fpvec[0], 0.996676642903, 1e-10);
    CHECK_CLOSE(Fpvec[1], 0.993739222842, 1e-10);
    CHECK_CLOSE(Fpvec[2], 0.990817393658, 1e-10);
    CHECK_CLOSE(Fpvec[3], 0.987911056819, 1e-10);
    CHECK_CLOSE(Fpvec[4], 0.985020114475, 1e-10);
    CHECK_CLOSE(Fpvec[5], 0.982144469445, 1e-10);
    CHECK_CLOSE(Fpvec[6], 0.97928402522, 1e-10);
    CHECK_CLOSE(Fpvec[7], 0.976438685953, 1e-10);
    CHECK_CLOSE(Fpvec[8], 0.973608356454, 1e-10);
    CHECK_CLOSE(Fpvec[9], 0.97079294219, 1e-10);


    //test for F0(x) for intermediate x:


    //test for F0(x) for large x:
    xvec = linspace(20,100,10);
    Fpvec = zeros(10);
    for(uint i = 0; i < 10; i++){
      for(uint n = 100; n > 0; n--){
               Fpvec[i] = boys.downwardRecursion(xvec[i],n);
    }
      }
    CHECK_CLOSE(Fpvec[0], 0.19816636483, 1e-8);
    CHECK_CLOSE(Fpvec[1], 0.164884382227, 1e-10);
    CHECK_CLOSE(Fpvec[2], 0.144187209502, 1e-10);
    CHECK_CLOSE(Fpvec[3], 0.12973033818, 1e-10);
    CHECK_CLOSE(Fpvec[4], 0.118899818928, 1e-10);
    CHECK_CLOSE(Fpvec[5], 0.110395710425, 1e-10);
    CHECK_CLOSE(Fpvec[6], 0.103489008863, 1e-10);
    CHECK_CLOSE(Fpvec[7], 0.0977350491129, 1e-10);
    CHECK_CLOSE(Fpvec[8], 0.0928451600494, 1e-10);
    CHECK_CLOSE(Fpvec[9], 0.0886226925453, 1e-10);

}


TEST(boysFunction)
{
    Boys boys;

    //test for Fn(0):
    CHECK_CLOSE(boys.boysFunction(0,0), 1.0, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,1), 0.333333333333, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,2), 0.2, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,3), 0.142857142857, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,4), 0.111111111111, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,5), 0.0909090909091, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,6), 0.0769230769231, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,7), 0.0666666666667, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,8), 0.0588235294118, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0,9), 0.0526315789474, 1e-10);

    //test for F0(x) for small x:
    CHECK_CLOSE(boys.boysFunction(0.01,0), 0.996676642903, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0188888888889,0), 0.993739222842, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0277777777778,0), 0.990817393658, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0366666666667,0), 0.987911056819, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0455555555556,0), 0.985020114475, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0544444444444,0), 0.982144469445, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0633333333333,0), 0.97928402522, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0722222222222,0), 0.976438685953, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.0811111111111,0), 0.973608356454, 1e-10);
    CHECK_CLOSE(boys.boysFunction(0.09,0), 0.97079294219, 1e-10);


    //test for F0(x) for intermediate x:
    //    CHECK_CLOSE(boys.boysFunction(0.1,0), 1.0, 1e-10);

    //test for F0(x) for large x:
    CHECK_CLOSE(boys.boysFunction(20.0,0), 0.19816636483, 1e-10);
    CHECK_CLOSE(boys.boysFunction(28.8888888889,0), 0.164884382227, 1e-10);
    CHECK_CLOSE(boys.boysFunction(37.7777777778,0), 0.144187209502, 1e-10);
    CHECK_CLOSE(boys.boysFunction(46.6666666667,0), 0.12973033818, 1e-10);
    CHECK_CLOSE(boys.boysFunction(55.5555555556,0), 0.118899818928, 1e-10);
    CHECK_CLOSE(boys.boysFunction(64.4444444444,0), 0.110395710425, 1e-10);
    CHECK_CLOSE(boys.boysFunction(73.3333333333,0), 0.103489008863, 1e-10);
    CHECK_CLOSE(boys.boysFunction(82.2222222222,0), 0.0977350491129, 1e-10);
    CHECK_CLOSE(boys.boysFunction(91.1111111111,0), 0.0928451600494, 1e-10);
    CHECK_CLOSE(boys.boysFunction(100.0,0), 0.0886226925453, 1e-10);
}

int main()
{
    return UnitTest::RunAllTests();
}
