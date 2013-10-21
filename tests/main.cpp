#include <unittest++/UnitTest++.h>
#include <src/integrator/integrator.h>


TEST(GTOIntegration) {
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

int main()
{
    return UnitTest::RunAllTests();
}
