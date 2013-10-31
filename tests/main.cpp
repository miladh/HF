#include <unittest++/UnitTest++.h>
#include <src/integrator/integrator.h>
#include <src/math/boys.h>


TEST(GTOoverlapIntegral)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);

//    integrator.addPrimitives(new PrimitiveGTO(0.2, 1.0));
//    integrator.addPrimitives(new PrimitiveGTO(0.3, 1.0));

    integrator.setMaxAngularMomentum(1);

    integrator.setupE();

    CHECK_CLOSE(1.191723635809e-01, integrator.overlapIntegral(0,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(2.764798835076e-01, integrator.overlapIntegral(0,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(7.606056933184e-01, integrator.overlapIntegral(0,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(4.290205088911e-02, integrator.overlapIntegral(0,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(9.953275806273e-02, integrator.overlapIntegral(0,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(1.346171019009e-01, integrator.overlapIntegral(0,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(1.191723635809e-01, integrator.overlapIntegral(0,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(2.764798835076e-01, integrator.overlapIntegral(0,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(4.290205088911e-02, integrator.overlapIntegral(0,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(2.383447271617e-01, integrator.overlapIntegral(0,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(0,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(-8.429776310255e-01, integrator.overlapIntegral(0,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(-2.093948045733e+00, integrator.overlapIntegral(0,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(0,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(-3.034719471692e-01, integrator.overlapIntegral(0,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(-4.684675146152e-01, integrator.overlapIntegral(0,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(0,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(-8.429776310255e-01, integrator.overlapIntegral(0,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(0,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(-8.294396505227e-01, integrator.overlapIntegral(0,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(1.562397355490e+00, integrator.overlapIntegral(0,0,2,0,0,0), 1e-5);
    CHECK_CLOSE(2.795322214215e+00, integrator.overlapIntegral(0,0,2,0,0,1), 1e-5);
    CHECK_CLOSE(6.361589630418e+00, integrator.overlapIntegral(0,0,2,0,0,2), 1e-5);
    CHECK_CLOSE(5.624630479765e-01, integrator.overlapIntegral(0,0,2,0,1,0), 1e-5);
    CHECK_CLOSE(1.006315997117e+00, integrator.overlapIntegral(0,0,2,0,1,1), 1e-5);
    CHECK_CLOSE(1.764884052762e+00, integrator.overlapIntegral(0,0,2,0,2,0), 1e-5);
    CHECK_CLOSE(1.562397355490e+00, integrator.overlapIntegral(0,0,2,1,0,0), 1e-5);
    CHECK_CLOSE(2.795322214215e+00, integrator.overlapIntegral(0,0,2,1,0,1), 1e-5);
    CHECK_CLOSE(5.624630479765e-01, integrator.overlapIntegral(0,0,2,1,1,0), 1e-5);
    CHECK_CLOSE(3.124794710981e+00, integrator.overlapIntegral(0,0,2,2,0,0), 1e-5);
    CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(0,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(0,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(-4.107270743920e-01, integrator.overlapIntegral(0,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(9.600525610073e-02, integrator.overlapIntegral(0,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(2.227321941537e-01, integrator.overlapIntegral(0,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(1.311086675171e-02, integrator.overlapIntegral(0,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(0,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(0,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(9.600525610073e-02, integrator.overlapIntegral(0,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(-1.287061526673e-01, integrator.overlapIntegral(0,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(0,1,1,0,0,0), 1e-5);
    CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(0,1,1,0,0,1), 1e-5);
    CHECK_CLOSE(1.130731944696e+00, integrator.overlapIntegral(0,1,1,0,0,2), 1e-5);
    CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(0,1,1,0,1,0), 1e-5);
    CHECK_CLOSE(-6.791027795542e-01, integrator.overlapIntegral(0,1,1,0,1,1), 1e-5);
    CHECK_CLOSE(-4.562581629595e-02, integrator.overlapIntegral(0,1,1,0,2,0), 1e-5);
    CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(0,1,1,1,0,0), 1e-5);
    CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(0,1,1,1,0,1), 1e-5);
    CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(0,1,1,1,1,0), 1e-5);
    CHECK_CLOSE(4.478974112823e-01, integrator.overlapIntegral(0,1,1,2,0,0), 1e-5);
    CHECK_CLOSE(1.539230248010e-01, integrator.overlapIntegral(0,2,0,0,0,0), 1e-5);
    CHECK_CLOSE(3.571014175384e-01, integrator.overlapIntegral(0,2,0,0,0,1), 1e-5);
    CHECK_CLOSE(9.823983134901e-01, integrator.overlapIntegral(0,2,0,0,0,2), 1e-5);
    CHECK_CLOSE(-7.329386373895e-02, integrator.overlapIntegral(0,2,0,0,1,0), 1e-5);
    CHECK_CLOSE(-1.700417638744e-01, integrator.overlapIntegral(0,2,0,0,1,1), 1e-5);
    CHECK_CLOSE(3.195477460565e-01, integrator.overlapIntegral(0,2,0,0,2,0), 1e-5);
    CHECK_CLOSE(1.539230248010e-01, integrator.overlapIntegral(0,2,0,1,0,0), 1e-5);
    CHECK_CLOSE(3.571014175384e-01, integrator.overlapIntegral(0,2,0,1,0,1), 1e-5);
    CHECK_CLOSE(-7.329386373895e-02, integrator.overlapIntegral(0,2,0,1,1,0), 1e-5);
    CHECK_CLOSE(3.078460496021e-01, integrator.overlapIntegral(0,2,0,2,0,0), 1e-5);
    CHECK_CLOSE(-1.787585453713e-01, integrator.overlapIntegral(1,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(-4.147198252614e-01, integrator.overlapIntegral(1,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(-1.140908539978e+00, integrator.overlapIntegral(1,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(-6.435307633366e-02, integrator.overlapIntegral(1,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-1.492991370941e-01, integrator.overlapIntegral(1,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-2.019256528514e-01, integrator.overlapIntegral(1,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(-5.958618179043e-02, integrator.overlapIntegral(1,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-1.382399417538e-01, integrator.overlapIntegral(1,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-2.145102544455e-02, integrator.overlapIntegral(1,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-1.191723635809e-01, integrator.overlapIntegral(1,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(6.220797378920e-01, integrator.overlapIntegral(1,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(1.264466446538e+00, integrator.overlapIntegral(1,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(3.140922068599e+00, integrator.overlapIntegral(1,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(1,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(4.552079207538e-01, integrator.overlapIntegral(1,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(7.027012719229e-01, integrator.overlapIntegral(1,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(2.073599126307e-01, integrator.overlapIntegral(1,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(4.214888155128e-01, integrator.overlapIntegral(1,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(7.464956854705e-02, integrator.overlapIntegral(1,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(4.147198252614e-01, integrator.overlapIntegral(1,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(9.652961450049e-02, integrator.overlapIntegral(1,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(2.239487056411e-01, integrator.overlapIntegral(1,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(6.160906115879e-01, integrator.overlapIntegral(1,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(-1.440078841511e-01, integrator.overlapIntegral(1,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(-3.340982912306e-01, integrator.overlapIntegral(1,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(-1.966630012757e-02, integrator.overlapIntegral(1,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(3.217653816683e-02, integrator.overlapIntegral(1,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(7.464956854705e-02, integrator.overlapIntegral(1,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(-4.800262805037e-02, integrator.overlapIntegral(1,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(6.435307633366e-02, integrator.overlapIntegral(1,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(3.873101816378e-01, integrator.overlapIntegral(2,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(8.985596213996e-01, integrator.overlapIntegral(2,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.471968503285e+00, integrator.overlapIntegral(2,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(1.394316653896e-01, integrator.overlapIntegral(2,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(3.234814637039e-01, integrator.overlapIntegral(2,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(4.375055811780e-01, integrator.overlapIntegral(2,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(2.979309089521e-02, integrator.overlapIntegral(2,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(6.911997087690e-02, integrator.overlapIntegral(2,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(1.072551272228e-02, integrator.overlapIntegral(2,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(2.979309089521e-01, integrator.overlapIntegral(2,0,0,2,0,0), 1e-5);
}
TEST(GTOkineticIntegral)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);

//    integrator.addPrimitives(new PrimitiveGTO(0.2, 1.0));
//    integrator.addPrimitives(new PrimitiveGTO(0.3, 1.0));

    integrator.setMaxAngularMomentum(2);

    integrator.setupE();

    CHECK_CLOSE(-9.678702680582e-02, integrator.kineticIntegral(0,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(-1.581907301477e-01, integrator.kineticIntegral(0,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(-3.289130997844e-01, integrator.kineticIntegral(0,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(-2.454683743671e-02, integrator.kineticIntegral(0,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-3.306080091812e-02, integrator.kineticIntegral(0,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-1.209847292592e-01, integrator.kineticIntegral(0,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(-6.818565954642e-02, integrator.kineticIntegral(0,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-9.183555810588e-02, integrator.kineticIntegral(0,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-1.425034522333e-02, integrator.kineticIntegral(0,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-1.554388972658e-01, integrator.kineticIntegral(0,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(2.372860952215e-01, integrator.kineticIntegral(0,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(2.514020826620e-01, integrator.kineticIntegral(0,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(3.256885959785e-01, integrator.kineticIntegral(0,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(4.959120137717e-02, integrator.kineticIntegral(0,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(1.767148243772e-02, integrator.kineticIntegral(0,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(3.085946543142e-01, integrator.kineticIntegral(0,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(1.377533371588e-01, integrator.kineticIntegral(0,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(4.908745121589e-02, integrator.kineticIntegral(0,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(0,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(3.418618463595e-01, integrator.kineticIntegral(0,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(-6.190706910076e-01, integrator.kineticIntegral(0,0,2,0,0,0), 1e-5);
    CHECK_CLOSE(-2.907944812827e-01, integrator.kineticIntegral(0,0,2,0,0,1), 1e-5);
    CHECK_CLOSE(1.109128468812e-01, integrator.kineticIntegral(0,0,2,0,0,2), 1e-5);
    CHECK_CLOSE(-8.787431724837e-02, integrator.kineticIntegral(0,0,2,0,1,0), 1e-5);
    CHECK_CLOSE(1.368298260464e-01, integrator.kineticIntegral(0,0,2,0,1,1), 1e-5);
    CHECK_CLOSE(-8.520922147503e-01, integrator.kineticIntegral(0,0,2,0,2,0), 1e-5);
    CHECK_CLOSE(-2.440953256899e-01, integrator.kineticIntegral(0,0,2,1,0,0), 1e-5);
    CHECK_CLOSE(3.800828501289e-01, integrator.kineticIntegral(0,0,2,1,0,1), 1e-5);
    CHECK_CLOSE(4.711681426600e-02, integrator.kineticIntegral(0,0,2,1,1,0), 1e-5);
    CHECK_CLOSE(-7.381742282583e-01, integrator.kineticIntegral(0,0,2,2,0,0), 1e-5);
    CHECK_CLOSE(3.682025615507e-02, integrator.kineticIntegral(0,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(4.959120137717e-02, integrator.kineticIntegral(0,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(7.903857602951e-02, integrator.kineticIntegral(0,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(-6.049047312582e-02, integrator.kineticIntegral(0,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(-8.688217105502e-02, integrator.kineticIntegral(0,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(1.938468694693e-02, integrator.kineticIntegral(0,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(2.137551783499e-02, integrator.kineticIntegral(0,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(0,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(-3.744921166165e-02, integrator.kineticIntegral(0,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(5.304752788336e-02, integrator.kineticIntegral(0,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(-7.438680206576e-02, integrator.kineticIntegral(0,1,1,0,0,0), 1e-5);
    CHECK_CLOSE(-2.650722365658e-02, integrator.kineticIntegral(0,1,1,0,0,1), 1e-5);
    CHECK_CLOSE(9.550382489858e-02, integrator.kineticIntegral(0,1,1,0,0,2), 1e-5);
    CHECK_CLOSE(1.303232565825e-01, integrator.kineticIntegral(0,1,1,0,1,0), 1e-5);
    CHECK_CLOSE(7.887481505265e-02, integrator.kineticIntegral(0,1,1,0,1,1), 1e-5);
    CHECK_CLOSE(-7.840890648636e-02, integrator.kineticIntegral(0,1,1,0,2,0), 1e-5);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(0,1,1,1,0,0), 1e-5);
    CHECK_CLOSE(8.274267732432e-02, integrator.kineticIntegral(0,1,1,1,0,1), 1e-5);
    CHECK_CLOSE(5.013966668721e-02, integrator.kineticIntegral(0,1,1,1,1,0), 1e-5);
    CHECK_CLOSE(-7.711001832636e-02, integrator.kineticIntegral(0,1,1,2,0,0), 1e-5);
    CHECK_CLOSE(-1.512318573258e-01, integrator.kineticIntegral(0,2,0,0,0,0), 1e-5);
    CHECK_CLOSE(-2.651535687867e-01, integrator.kineticIntegral(0,2,0,0,0,1), 1e-5);
    CHECK_CLOSE(-5.921817515938e-01, integrator.kineticIntegral(0,2,0,0,0,2), 1e-5);
    CHECK_CLOSE(1.606516375488e-03, integrator.kineticIntegral(0,2,0,0,1,0), 1e-5);
    CHECK_CLOSE(-3.708290533872e-02, integrator.kineticIntegral(0,2,0,0,1,1), 1e-5);
    CHECK_CLOSE(-2.565122085492e-01, integrator.kineticIntegral(0,2,0,0,2,0), 1e-5);
    CHECK_CLOSE(-1.142903313736e-01, integrator.kineticIntegral(0,2,0,1,0,0), 1e-5);
    CHECK_CLOSE(-1.794492285775e-01, integrator.kineticIntegral(0,2,0,1,0,1), 1e-5);
    CHECK_CLOSE(-1.598401092187e-02, integrator.kineticIntegral(0,2,0,1,1,0), 1e-5);
    CHECK_CLOSE(-2.532083467153e-01, integrator.kineticIntegral(0,2,0,2,0,0), 1e-5);
    CHECK_CLOSE(1.022784893196e-01, integrator.kineticIntegral(1,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(1.377533371588e-01, integrator.kineticIntegral(1,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.195516000820e-01, integrator.kineticIntegral(1,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(2.137551783499e-02, integrator.kineticIntegral(1,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(1.375940847459e-02, integrator.kineticIntegral(1,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(1.330149372044e-01, integrator.kineticIntegral(1,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(-8.809221115897e-03, integrator.kineticIntegral(1,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-5.361497900979e-02, integrator.kineticIntegral(1,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-8.319565708416e-03, integrator.kineticIntegral(1,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(6.818565954642e-02, integrator.kineticIntegral(1,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(-2.066300057382e-01, integrator.kineticIntegral(1,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(-7.363117682384e-02, integrator.kineticIntegral(1,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(2.652884024961e-01, integrator.kineticIntegral(1,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(1,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(8.274267732432e-02, integrator.kineticIntegral(1,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(-2.942436762098e-01, integrator.kineticIntegral(1,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(8.042246851468e-02, integrator.kineticIntegral(1,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(2.789282215612e-01, integrator.kineticIntegral(1,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(4.686798511658e-02, integrator.kineticIntegral(1,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(-1.377533371588e-01, integrator.kineticIntegral(1,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(-3.206327675248e-02, integrator.kineticIntegral(1,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(-2.063911271189e-02, integrator.kineticIntegral(1,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.930388273683e-02, integrator.kineticIntegral(1,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(5.617381749247e-02, integrator.kineticIntegral(1,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(5.013966668720e-02, integrator.kineticIntegral(1,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(-3.379694245101e-02, integrator.kineticIntegral(1,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(1.247934856262e-02, integrator.kineticIntegral(1,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(4.686798511658e-02, integrator.kineticIntegral(1,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(-1.583728636544e-02, integrator.kineticIntegral(1,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(-2.137551783499e-02, integrator.kineticIntegral(1,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(-2.287537353407e-01, integrator.kineticIntegral(2,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(-3.150543568546e-01, integrator.kineticIntegral(2,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(-5.213314751101e-01, integrator.kineticIntegral(2,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(-4.888774502915e-02, integrator.kineticIntegral(2,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-3.578401717871e-02, integrator.kineticIntegral(2,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-2.962760567236e-01, integrator.kineticIntegral(2,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(-1.704641488661e-02, integrator.kineticIntegral(2,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-2.295888952647e-02, integrator.kineticIntegral(2,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-3.562586305832e-03, integrator.kineticIntegral(2,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-3.468392469657e-01, integrator.kineticIntegral(2,0,0,2,0,0), 1e-5);
}
TEST(GTOnuclearAttractionIntegral)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec posC = {2.3,0.9,3.2};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setCorePositionC(posC);

//    integrator.addPrimitives(new PrimitiveGTO(0.2, 1.0));
//    integrator.addPrimitives(new PrimitiveGTO(0.3, 1.0));

    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);
    integrator.setMaxAngularMomentum(2);

    integrator.setupE();


    CHECK_CLOSE(2.788948987251e-02, integrator.nuclearAttractionIntegral(0,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(6.971203468743e-02, integrator.nuclearAttractionIntegral(0,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.024071525839e-01, integrator.nuclearAttractionIntegral(0,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(8.727033700014e-03, integrator.nuclearAttractionIntegral(0,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(2.134361291529e-02, integrator.nuclearAttractionIntegral(0,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(2.921666495443e-02, integrator.nuclearAttractionIntegral(0,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(3.185957751329e-02, integrator.nuclearAttractionIntegral(0,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(8.105746642202e-02, integrator.nuclearAttractionIntegral(0,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(9.596523510045e-03, integrator.nuclearAttractionIntegral(0,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(6.388444040338e-02, integrator.nuclearAttractionIntegral(0,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(-9.204700718547e-02, integrator.nuclearAttractionIntegral(0,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(-2.019226499202e-01, integrator.nuclearAttractionIntegral(0,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(-5.276399683274e-01, integrator.nuclearAttractionIntegral(0,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(-2.927318225377e-02, integrator.nuclearAttractionIntegral(0,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(-6.299787002318e-02, integrator.nuclearAttractionIntegral(0,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(-9.718105595370e-02, integrator.nuclearAttractionIntegral(0,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(-1.037280861539e-01, integrator.nuclearAttractionIntegral(0,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(-2.312309453843e-01, integrator.nuclearAttractionIntegral(0,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(-3.202910466576e-02, integrator.nuclearAttractionIntegral(0,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(-2.073449397904e-01, integrator.nuclearAttractionIntegral(0,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(3.319499900436e-01, integrator.nuclearAttractionIntegral(0,0,2,0,0,0), 1e-5);
    CHECK_CLOSE(6.435114042344e-01, integrator.nuclearAttractionIntegral(0,0,2,0,0,1), 1e-5);
    CHECK_CLOSE(1.536931448007e+00, integrator.nuclearAttractionIntegral(0,0,2,0,0,2), 1e-5);
    CHECK_CLOSE(1.067865861209e-01, integrator.nuclearAttractionIntegral(0,0,2,0,1,0), 1e-5);
    CHECK_CLOSE(2.033153544029e-01, integrator.nuclearAttractionIntegral(0,0,2,0,1,1), 1e-5);
    CHECK_CLOSE(3.524622701603e-01, integrator.nuclearAttractionIntegral(0,0,2,0,2,0), 1e-5);
    CHECK_CLOSE(3.703919381580e-01, integrator.nuclearAttractionIntegral(0,0,2,1,0,0), 1e-5);
    CHECK_CLOSE(7.292169308884e-01, integrator.nuclearAttractionIntegral(0,0,2,1,0,1), 1e-5);
    CHECK_CLOSE(1.162963233448e-01, integrator.nuclearAttractionIntegral(0,0,2,1,1,0), 1e-5);
    CHECK_CLOSE(7.390872806284e-01, integrator.nuclearAttractionIntegral(0,0,2,2,0,0), 1e-5);
    CHECK_CLOSE(-1.637350724302e-02, integrator.nuclearAttractionIntegral(0,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(-4.139721853567e-02, integrator.nuclearAttractionIntegral(0,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(-1.213713540367e-01, integrator.nuclearAttractionIntegral(0,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(2.136233458775e-02, integrator.nuclearAttractionIntegral(0,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(5.306634838230e-02, integrator.nuclearAttractionIntegral(0,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(-1.697963144320e-04, integrator.nuclearAttractionIntegral(0,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(-1.907709578263e-02, integrator.nuclearAttractionIntegral(0,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(-4.932098923684e-02, integrator.nuclearAttractionIntegral(0,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(2.414126847830e-02, integrator.nuclearAttractionIntegral(0,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(-3.842342257899e-02, integrator.nuclearAttractionIntegral(0,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(5.356912442721e-02, integrator.nuclearAttractionIntegral(0,1,1,0,0,0), 1e-5);
    CHECK_CLOSE(1.187325132374e-01, integrator.nuclearAttractionIntegral(0,1,1,0,0,1), 1e-5);
    CHECK_CLOSE(3.128036711345e-01, integrator.nuclearAttractionIntegral(0,1,1,0,0,2), 1e-5);
    CHECK_CLOSE(-7.083519267815e-02, integrator.nuclearAttractionIntegral(0,1,1,0,1,0), 1e-5);
    CHECK_CLOSE(-1.544897767601e-01, integrator.nuclearAttractionIntegral(0,1,1,0,1,1), 1e-5);
    CHECK_CLOSE(-3.894393296797e-04, integrator.nuclearAttractionIntegral(0,1,1,0,2,0), 1e-5);
    CHECK_CLOSE(6.132616876012e-02, integrator.nuclearAttractionIntegral(0,1,1,1,0,0), 1e-5);
    CHECK_CLOSE(1.386353605834e-01, integrator.nuclearAttractionIntegral(0,1,1,1,0,1), 1e-5);
    CHECK_CLOSE(-7.911527548440e-02, integrator.nuclearAttractionIntegral(0,1,1,1,1,0), 1e-5);
    CHECK_CLOSE(1.229240947242e-01, integrator.nuclearAttractionIntegral(0,1,1,2,0,0), 1e-5);
    CHECK_CLOSE(3.609849112824e-02, integrator.nuclearAttractionIntegral(0,2,0,0,0,0), 1e-5);
    CHECK_CLOSE(9.032384498820e-02, integrator.nuclearAttractionIntegral(0,2,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.625292648498e-01, integrator.nuclearAttractionIntegral(0,2,0,0,0,2), 1e-5);
    CHECK_CLOSE(-1.939589748931e-02, integrator.nuclearAttractionIntegral(0,2,0,0,1,0), 1e-5);
    CHECK_CLOSE(-4.913397190183e-02, integrator.nuclearAttractionIntegral(0,2,0,0,1,1), 1e-5);
    CHECK_CLOSE(6.878262296370e-02, integrator.nuclearAttractionIntegral(0,2,0,0,2,0), 1e-5);
    CHECK_CLOSE(4.131065513841e-02, integrator.nuclearAttractionIntegral(0,2,0,1,0,0), 1e-5);
    CHECK_CLOSE(1.052929737663e-01, integrator.nuclearAttractionIntegral(0,2,0,1,0,1), 1e-5);
    CHECK_CLOSE(-2.267402937768e-02, integrator.nuclearAttractionIntegral(0,2,0,1,1,0), 1e-5);
    CHECK_CLOSE(8.289710831960e-02, integrator.nuclearAttractionIntegral(0,2,0,2,0,0), 1e-5);
    CHECK_CLOSE(-3.786414780578e-02, integrator.nuclearAttractionIntegral(1,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(-9.322262110550e-02, integrator.nuclearAttractionIntegral(1,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(-2.671155215998e-01, integrator.nuclearAttractionIntegral(1,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(-1.222106053447e-02, integrator.nuclearAttractionIntegral(1,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-2.972830178046e-02, integrator.nuclearAttractionIntegral(1,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-4.026352276293e-02, integrator.nuclearAttractionIntegral(1,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(-1.576450345257e-02, integrator.nuclearAttractionIntegral(1,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(-3.945885129414e-02, integrator.nuclearAttractionIntegral(1,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(-4.918734877201e-03, integrator.nuclearAttractionIntegral(1,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-2.459437143524e-02, integrator.nuclearAttractionIntegral(1,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(1.263894353489e-01, integrator.nuclearAttractionIntegral(1,0,1,0,0,0), 1e-5);
    CHECK_CLOSE(2.735756798558e-01, integrator.nuclearAttractionIntegral(1,0,1,0,0,1), 1e-5);
    CHECK_CLOSE(7.071773603054e-01, integrator.nuclearAttractionIntegral(1,0,1,0,0,2), 1e-5);
    CHECK_CLOSE(4.115384967585e-02, integrator.nuclearAttractionIntegral(1,0,1,0,1,0), 1e-5);
    CHECK_CLOSE(8.802219023191e-02, integrator.nuclearAttractionIntegral(1,0,1,0,1,1), 1e-5);
    CHECK_CLOSE(1.350111738398e-01, integrator.nuclearAttractionIntegral(1,0,1,0,2,0), 1e-5);
    CHECK_CLOSE(5.197526800952e-02, integrator.nuclearAttractionIntegral(1,0,1,1,0,0), 1e-5);
    CHECK_CLOSE(1.145639876363e-01, integrator.nuclearAttractionIntegral(1,0,1,1,0,1), 1e-5);
    CHECK_CLOSE(1.638641395940e-02, integrator.nuclearAttractionIntegral(1,0,1,1,1,0), 1e-5);
    CHECK_CLOSE(8.278875254192e-02, integrator.nuclearAttractionIntegral(1,0,1,2,0,0), 1e-5);
    CHECK_CLOSE(2.185667243163e-02, integrator.nuclearAttractionIntegral(1,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(5.417205698627e-02, integrator.nuclearAttractionIntegral(1,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(1.560020091608e-01, integrator.nuclearAttractionIntegral(1,1,0,0,0,2), 1e-5);
    CHECK_CLOSE(-2.926456829930e-02, integrator.nuclearAttractionIntegral(1,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(-7.176178735649e-02, integrator.nuclearAttractionIntegral(1,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(-5.223967979758e-04, integrator.nuclearAttractionIntegral(1,1,0,0,2,0), 1e-5);
    CHECK_CLOSE(9.269318129877e-03, integrator.nuclearAttractionIntegral(1,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(2.337071697343e-02, integrator.nuclearAttractionIntegral(1,1,0,1,0,1), 1e-5);
    CHECK_CLOSE(-1.203714316117e-02, integrator.nuclearAttractionIntegral(1,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(1.401501778682e-02, integrator.nuclearAttractionIntegral(1,1,0,2,0,0), 1e-5);
    CHECK_CLOSE(7.889586550718e-02, integrator.nuclearAttractionIntegral(2,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(1.935977010010e-01, integrator.nuclearAttractionIntegral(2,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(5.534914541236e-01, integrator.nuclearAttractionIntegral(2,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(2.563391673303e-02, integrator.nuclearAttractionIntegral(2,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(6.217850538435e-02, integrator.nuclearAttractionIntegral(2,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(8.419480232293e-02, integrator.nuclearAttractionIntegral(2,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(1.481688684288e-02, integrator.nuclearAttractionIntegral(2,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(3.878852644576e-02, integrator.nuclearAttractionIntegral(2,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(4.176920693786e-03, integrator.nuclearAttractionIntegral(2,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(6.422210627967e-02, integrator.nuclearAttractionIntegral(2,0,0,2,0,0), 1e-5);

}
TEST(BoysfactorialFunctions)
{
    Boys boys(0);
    CHECK_EQUAL(boys.factorial(0),1);
    CHECK_EQUAL(boys.factorial(1),1);
    CHECK_EQUAL(boys.factorial(5),120);
    CHECK_EQUAL(boys.factorial(20),2432902008176640000);
    CHECK_EQUAL(boys.factorial(-1),1);

    CHECK_EQUAL(boys.doubleFactorial(0),1);
    CHECK_EQUAL(boys.doubleFactorial(1),1);
    CHECK_EQUAL(boys.doubleFactorial(5),15);
    CHECK_EQUAL(boys.doubleFactorial(-1),1);
}
TEST(boysFunction)
{
    //test for F0(x) for small x <= 50:
    //Recursion relation is not used!
    Boys boysF0_small(0);
    mat F0; rowvec xvec;
    xvec = linspace<rowvec>(0.01,50,20);
    F0 = zeros(20,1);
    for(uint i = 0; i < 20; i++){
        boysF0_small.evaluateBoysFunctions(xvec[i]);
        F0.row(i) =  boysF0_small.getBoysFunctions();
    }

    CHECK_CLOSE(F0(0,0), 9.9667664290336333e-01, 1e-10);
    CHECK_CLOSE(F0(1,0), 5.3357683580906246e-01, 1e-10);
    CHECK_CLOSE(F0(2,0), 3.8551956882369670e-01, 1e-10);
    CHECK_CLOSE(F0(3,0), 3.1522027004511449e-01, 1e-10);
    CHECK_CLOSE(F0(4,0), 2.7304989839411259e-01, 1e-10);
    CHECK_CLOSE(F0(5,0), 2.4424745291117503e-01, 1e-10);
    CHECK_CLOSE(F0(6,0), 2.2298057384240719e-01, 1e-10);
    CHECK_CLOSE(F0(7,0), 2.0644923632081960e-01, 1e-10);
    CHECK_CLOSE(F0(8,0), 1.9312212797082015e-01, 1e-10);
    CHECK_CLOSE(F0(9,0), 1.8208209208151274e-01, 1e-10);
    CHECK_CLOSE(F0(10,0), 1.7274188563417292e-01, 1e-10);
    CHECK_CLOSE(F0(11,0), 1.6470576997839928e-01, 1e-10);
    CHECK_CLOSE(F0(12,0), 1.5769603853464209e-01, 1e-10);
    CHECK_CLOSE(F0(13,0), 1.5151129820345471e-01, 1e-10);
    CHECK_CLOSE(F0(14,0), 1.4600146419605400e-01, 1e-10);
    CHECK_CLOSE(F0(15,0), 1.4105209097756166e-01, 1e-10);
    CHECK_CLOSE(F0(16,0), 1.3657418098484136e-01, 1e-10);
    CHECK_CLOSE(F0(17,0), 1.3249734289737117e-01, 1e-10);
    CHECK_CLOSE(F0(18,0), 1.2876507161025572e-01, 1e-10);
    CHECK_CLOSE(F0(19,0), 1.2533141373155002e-01, 1e-10);




    //test for F0(x) for large x > 50:
    //Recursion relation is not used!
    Boys boysF0_large(0);
    xvec = linspace<rowvec>(51,100,20);
    F0 = zeros(20,1);
    for(uint i = 0; i < 20; i++){
        boysF0_large.evaluateBoysFunctions(xvec[i]);
        F0.row(i) = boysF0_large.getBoysFunctions();
    }

    CHECK_CLOSE(F0(0,0), 1.2409659136408727e-01, 1e-10);
    CHECK_CLOSE(F0(1,0), 1.2107315290425182e-01, 1e-10);
    CHECK_CLOSE(F0(2,0), 1.1826045114934411e-01, 1e-10);
    CHECK_CLOSE(F0(3,0), 1.1563509029711123e-01, 1e-10);
    CHECK_CLOSE(F0(4,0), 1.1317715652582763e-01, 1e-10);
    CHECK_CLOSE(F0(5,0), 1.1086957884293665e-01, 1e-10);
    CHECK_CLOSE(F0(6,0), 1.0869762771661817e-01, 1e-10);
    CHECK_CLOSE(F0(7,0), 1.0664851770975842e-01, 1e-10);
    CHECK_CLOSE(F0(8,0), 1.0471108954515591e-01, 1e-10);
    CHECK_CLOSE(F0(9,0), 1.0287555349437555e-01, 1e-10);
    CHECK_CLOSE(F0(10,0), 1.0113328058446551e-01, 1e-10);
    CHECK_CLOSE(F0(11,0), 9.9476631436519997e-02, 1e-10);
    CHECK_CLOSE(F0(12,0), 9.7898814974335960e-02, 1e-10);
    CHECK_CLOSE(F0(13,0), 9.6393771031862641e-02, 1e-10);
    CHECK_CLOSE(F0(14,0), 9.4956072224460661e-02, 1e-10);
    CHECK_CLOSE(F0(15,0), 9.3580841456184838e-02, 1e-10);
    CHECK_CLOSE(F0(16,0), 9.2263682201430275e-02, 1e-10);
    CHECK_CLOSE(F0(17,0), 9.1000619287057175e-02, 1e-10);
    CHECK_CLOSE(F0(18,0), 8.9788048355705335e-02, 1e-10);
    CHECK_CLOSE(F0(19,0), 8.8622692545275800e-02, 1e-10);
}
TEST(BoysDownwardrecursionFunction)
{

    //test for Fn(0):
    Boys boysFn_0(20);
    boysFn_0.evaluateBoysFunctions(0);
    rowvec Fn = boysFn_0.getBoysFunctions();
    CHECK_CLOSE(Fn[0], 1.0000000000000000e+00, 1e-10);
    CHECK_CLOSE(Fn[1], 3.3333333333333331e-01, 1e-10);
    CHECK_CLOSE(Fn[2], 2.0000000000000001e-01, 1e-10);
    CHECK_CLOSE(Fn[3], 1.4285714285714285e-01, 1e-10);
    CHECK_CLOSE(Fn[4], 1.1111111111111110e-01, 1e-10);
    CHECK_CLOSE(Fn[5], 9.0909090909090912e-02, 1e-10);
    CHECK_CLOSE(Fn[6], 7.6923076923076927e-02, 1e-10);
    CHECK_CLOSE(Fn[7], 6.6666666666666666e-02, 1e-10);
    CHECK_CLOSE(Fn[8], 5.8823529411764705e-02, 1e-10);
    CHECK_CLOSE(Fn[9], 5.2631578947368418e-02, 1e-10);
    CHECK_CLOSE(Fn[10], 4.7619047619047616e-02, 1e-10);
    CHECK_CLOSE(Fn[11], 4.3478260869565216e-02, 1e-10);
    CHECK_CLOSE(Fn[12], 4.0000000000000001e-02, 1e-10);
    CHECK_CLOSE(Fn[13], 3.7037037037037035e-02, 1e-10);
    CHECK_CLOSE(Fn[14], 3.4482758620689655e-02, 1e-10);
    CHECK_CLOSE(Fn[15], 3.2258064516129031e-02, 1e-10);
    CHECK_CLOSE(Fn[16], 3.0303030303030304e-02, 1e-10);
    CHECK_CLOSE(Fn[17], 2.8571428571428571e-02, 1e-10);
    CHECK_CLOSE(Fn[18], 2.7027027027027029e-02, 1e-10);
    CHECK_CLOSE(Fn[19], 2.5641025641025640e-02, 1e-10);
    CHECK_CLOSE(Fn[20], 2.4390243902439025e-02, 1e-9);




    //test for F0(x) for small x <= 50:
    // n goes from nMax to 0
    Boys boysF0_small(15);
    mat F0; rowvec xvec;
    xvec = linspace<rowvec>(0.01,50,20);
    F0 = zeros(20,16);
    for(uint i = 0; i < 20; i++){
        boysF0_small.evaluateBoysFunctions(xvec[i]);
        F0.row(i) =  boysF0_small.getBoysFunctions();
    }

    CHECK_CLOSE(F0(0,0), 9.9667664290336333e-01, 1e-10);
    CHECK_CLOSE(F0(1,0), 5.3357683580906246e-01, 1e-10);
    CHECK_CLOSE(F0(2,0), 3.8551956882369670e-01, 1e-10);
    CHECK_CLOSE(F0(3,0), 3.1522027004511449e-01, 1e-10);
    CHECK_CLOSE(F0(4,0), 2.7304989839411259e-01, 1e-10);
    CHECK_CLOSE(F0(5,0), 2.4424745291117503e-01, 1e-10);
    CHECK_CLOSE(F0(6,0), 2.2298057384240719e-01, 1e-10);
    CHECK_CLOSE(F0(7,0), 2.0644923632081960e-01, 1e-10);
    CHECK_CLOSE(F0(8,0), 1.9312212797082015e-01, 1e-10);
    CHECK_CLOSE(F0(9,0), 1.8208209208151274e-01, 1e-10);
    CHECK_CLOSE(F0(10,0), 1.7274188563417292e-01, 1e-10);
    CHECK_CLOSE(F0(11,0), 1.6470576997839928e-01, 1e-10);
    CHECK_CLOSE(F0(12,0), 1.5769603853464209e-01, 1e-10);
    CHECK_CLOSE(F0(13,0), 1.5151129820345471e-01, 1e-10);
    CHECK_CLOSE(F0(14,0), 1.4600146419605400e-01, 1e-10);
    CHECK_CLOSE(F0(15,0), 1.4105209097756166e-01, 1e-10);
    CHECK_CLOSE(F0(16,0), 1.3657418098484136e-01, 1e-10);
    CHECK_CLOSE(F0(17,0), 1.3249734289737117e-01, 1e-10);
    CHECK_CLOSE(F0(18,0), 1.2876507161025572e-01, 1e-10);
    CHECK_CLOSE(F0(19,0), 1.2533141373155002e-01, 1e-10);





    //test for F0(x) for large x > 50:
    // n goes from nMax to 0
    Boys boysF0_large(15);
    xvec = linspace<rowvec>(51,100,20);
    F0 = zeros(20,16);
    for(uint i = 0; i < 20; i++){
     boysF0_large.evaluateBoysFunctions(xvec[i]);
     F0.row(i) = boysF0_large.getBoysFunctions();
    }

    CHECK_CLOSE(F0(0,0), 1.2409659136408727e-01, 1e-9); //OBS!!
    CHECK_CLOSE(F0(1,0), 1.2107315290425182e-01, 1e-10);
    CHECK_CLOSE(F0(2,0), 1.1826045114934411e-01, 1e-10);
    CHECK_CLOSE(F0(3,0), 1.1563509029711123e-01, 1e-10);
    CHECK_CLOSE(F0(4,0), 1.1317715652582763e-01, 1e-10);
    CHECK_CLOSE(F0(5,0), 1.1086957884293665e-01, 1e-10);
    CHECK_CLOSE(F0(6,0), 1.0869762771661817e-01, 1e-10);
    CHECK_CLOSE(F0(7,0), 1.0664851770975842e-01, 1e-10);
    CHECK_CLOSE(F0(8,0), 1.0471108954515591e-01, 1e-10);
    CHECK_CLOSE(F0(9,0), 1.0287555349437555e-01, 1e-10);
    CHECK_CLOSE(F0(10,0), 1.0113328058446551e-01, 1e-10);
    CHECK_CLOSE(F0(11,0), 9.9476631436519997e-02, 1e-10);
    CHECK_CLOSE(F0(12,0), 9.7898814974335960e-02, 1e-10);
    CHECK_CLOSE(F0(13,0), 9.6393771031862641e-02, 1e-10);
    CHECK_CLOSE(F0(14,0), 9.4956072224460661e-02, 1e-10);
    CHECK_CLOSE(F0(15,0), 9.3580841456184838e-02, 1e-10);
    CHECK_CLOSE(F0(16,0), 9.2263682201430275e-02, 1e-10);
    CHECK_CLOSE(F0(17,0), 9.1000619287057175e-02, 1e-10);
    CHECK_CLOSE(F0(18,0), 8.9788048355705335e-02, 1e-10);
    CHECK_CLOSE(F0(19,0), 8.8622692545275800e-02, 1e-10);
}

int main()
{
    return UnitTest::RunAllTests();
}
