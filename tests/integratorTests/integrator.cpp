#include <unittest++/UnitTest++.h>
#include <integrator/integrator.h>

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

TEST(GTOoverlapIntegral)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);

    integrator.setMaxAngularMomentum(2);

    integrator.updateHermiteCoefficients(true, false);

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

TEST(GTOoverlapIntegral_derivative)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec Sab  = {0, 0, 0};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);

    integrator.setMaxAngularMomentum(1);

    integrator.updateHermiteCoefficients(true, false);
    integrator.updateHermiteCoefficients_derivative(true,false);

    Sab  = {-7.150341814851e-02,-2.574123053346e-02,-1.658879301045e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,0,0),3, 1e-5);

    Sab  = {-1.658879301045e-01,-5.971965483764e-02,-3.371910524102e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,0,1),3, 1e-5);

    Sab  = {-4.563634159911e-01,-1.642908297568e-01,-8.375792182932e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,0,2),3, 1e-5);

    Sab  = {-2.574123053346e-02,3.840210244029e-02,-5.971965483764e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,1,0),3, 1e-5);

    Sab  = {-5.971965483764e-02,8.909287766148e-02,-1.213887788677e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,1,1),3, 1e-5);

    Sab  = {-8.077026114056e-02,5.244346700684e-03,-1.873870058461e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,0,2,0),3, 1e-5);

    Sab  = {-2.383447271617e-02,-2.574123053346e-02,-1.658879301045e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,1,0,0),3, 1e-5);

    Sab  = {-5.529597670152e-02,-5.971965483764e-02,-3.371910524102e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,1,0,1),3, 1e-5);

    Sab  = {-8.580410177821e-03,3.840210244029e-02,-5.971965483764e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,1,1,0),3, 1e-5);

    Sab  = {-4.766894543234e-02,-5.148246106693e-02,-3.317758602091e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,0,2,0,0),3, 1e-5);

    Sab  = {2.488318951568e-01,8.957948225645e-02,5.057865786153e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,0,0),3, 1e-5);

    Sab  = {5.057865786153e-01,1.820831683015e-01,8.416490021784e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,0,1),3, 1e-5);

    Sab  = {1.256368827440e+00,4.522927778783e-01,1.784030158849e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,0,2),3, 1e-5);

    Sab  = {8.957948225646e-02,-1.336393164922e-01,1.820831683015e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,1,0),3, 1e-5);

    Sab  = {1.820831683015e-01,-2.716411118217e-01,3.029936407842e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,1,1),3, 1e-5);

    Sab  = {2.810805087691e-01,-1.825032651838e-02,5.713365192039e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,0,2,0),3, 1e-5);

    Sab  = {8.294396505227e-02,8.957948225646e-02,5.057865786153e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,1,0,0),3, 1e-5);

    Sab  = {1.685955262051e-01,1.820831683015e-01,8.416490021784e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,1,0,1),3, 1e-5);

    Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,1,1,0),3, 1e-5);

    Sab  = {1.658879301045e-01,1.791589645129e-01,1.011573157231e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,1,2,0,0),3, 1e-5);

    Sab  = {-9.374384132942e-01,-3.374778287859e-01,-1.677193328529e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,0,0),3, 1e-5);

    Sab  = {-1.677193328529e+00,-6.037895982704e-01,-2.254556422760e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,0,1),3, 1e-5);

    Sab  = {-3.816953778251e+00,-1.374103360170e+00,-4.106337339291e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,0,2),3, 1e-5);

    Sab  = {-3.374778287859e-01,5.034669238332e-01,-6.037895982704e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,1,0),3, 1e-5);

    Sab  = {-6.037895982704e-01,9.007646303086e-01,-8.116403121938e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,1,1),3, 1e-5);

    Sab  = {-1.058930431657e+00,6.875548298465e-02,-1.894557583906e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,0,2,0),3, 1e-5);

    Sab  = {-3.124794710981e-01,-3.374778287859e-01,-1.677193328529e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,1,0,0),3, 1e-5);

    Sab  = {-5.590644428430e-01,-6.037895982704e-01,-2.254556422760e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,1,0,1),3, 1e-5);

    Sab  = {-1.124926095953e-01,5.034669238332e-01,-6.037895982704e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,1,1,0),3, 1e-5);

    Sab  = {-6.249589421962e-01,-6.749556575719e-01,-3.354386657058e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,0,2,2,0,0),3, 1e-5);

    Sab  = {3.861184580020e-02,-5.760315366044e-02,8.957948225646e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,0,0),3, 1e-5);

    Sab  = {8.957948225646e-02,-1.336393164922e-01,1.820831683015e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,0,1),3, 1e-5);

    Sab  = {2.464362446352e-01,-3.676463679224e-01,4.522927778783e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,0,2),3, 1e-5);

    Sab  = {-5.760315366044e-02,-7.221959638469e-02,-1.336393164922e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,1,0),3, 1e-5);

    Sab  = {-1.336393164922e-01,-1.675494636125e-01,-2.716411118217e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,1,1),3, 1e-5);

    Sab  = {-7.866520051026e-03,-6.798003478342e-03,-1.825032651838e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,0,2,0),3, 1e-5);

    Sab  = {1.287061526673e-02,-5.760315366044e-02,8.957948225646e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,1,0,0),3, 1e-5);

    Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,1,0,1),3, 1e-5);

    Sab  = {-1.920105122015e-02,-7.221959638469e-02,-1.336393164922e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,1,1,0),3, 1e-5);

    Sab  = {2.574123053346e-02,-1.152063073209e-01,1.791589645129e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,0,2,0,0),3, 1e-5);

    Sab  = {-1.343692233847e-01,2.004589747383e-01,-2.731247524523e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,0,0),3, 1e-5);

    Sab  = {-2.731247524523e-01,4.074616677325e-01,-4.544904611764e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,0,1),3, 1e-5);

    Sab  = {-6.784391668175e-01,1.012130727385e+00,-9.633762857784e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,0,2),3, 1e-5);

    Sab  = {2.004589747383e-01,2.513241954187e-01,4.074616677325e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,1,0),3, 1e-5);

    Sab  = {4.074616677325e-01,5.108525369867e-01,6.780324361549e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,1,1),3, 1e-5);

    Sab  = {2.737548977757e-02,2.365705210462e-02,5.564461623294e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,0,2,0),3, 1e-5);

    Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,1,0,0),3, 1e-5);

    Sab  = {-9.104158415075e-02,4.074616677325e-01,-4.544904611764e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,1,0,1),3, 1e-5);

    Sab  = {6.681965824611e-02,2.513241954187e-01,4.074616677325e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,1,1,0),3, 1e-5);

    Sab  = {-8.957948225645e-02,4.009179494767e-01,-5.462495049045e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,1,1,2,0,0),3, 1e-5);

    Sab  = {-9.235381488062e-02,4.397631824337e-02,-2.142608505230e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,0,0),3, 1e-5);

    Sab  = {-2.142608505230e-01,1.020250583246e-01,-4.355159632930e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,0,1),3, 1e-5);

    Sab  = {-5.894389880941e-01,2.806744535565e-01,-1.081817318347e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,0,2),3, 1e-5);

    Sab  = {4.397631824337e-02,-3.780562283286e-02,1.020250583246e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,1,0),3, 1e-5);

    Sab  = {1.020250583246e-01,-8.770904497223e-02,2.073805898175e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,1,1),3, 1e-5);

    Sab  = {-1.917286476339e-01,-1.433904442414e-01,-4.448104625106e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,0,2,0),3, 1e-5);

    Sab  = {-3.078460496021e-02,4.397631824337e-02,-2.142608505230e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,1,0,0),3, 1e-5);

    Sab  = {-7.142028350768e-02,1.020250583246e-01,-4.355159632930e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,1,0,1),3, 1e-5);

    Sab  = {1.465877274779e-02,-3.780562283286e-02,1.020250583246e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,1,1,0),3, 1e-5);

    Sab  = {-6.156920992041e-02,8.795263648674e-02,-4.285217010461e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(0,2,0,2,0,0),3, 1e-5);

    Sab  = {3.575170907426e-02,3.861184580020e-02,2.488318951568e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,0,0),3, 1e-5);

    Sab  = {8.294396505227e-02,8.957948225646e-02,5.057865786153e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,0,1),3, 1e-5);

    Sab  = {2.281817079955e-01,2.464362446352e-01,1.256368827440e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,0,2),3, 1e-5);

    Sab  = {1.287061526673e-02,-5.760315366044e-02,8.957948225646e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,1,0),3, 1e-5);

    Sab  = {2.985982741882e-02,-1.336393164922e-01,1.820831683015e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,1,1),3, 1e-5);

    Sab  = {4.038513057028e-02,-7.866520051026e-03,2.810805087691e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,0,2,0),3, 1e-5);

    Sab  = {-1.072551272228e-01,1.287061526673e-02,8.294396505227e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,1,0,0),3, 1e-5);

    Sab  = {-2.488318951568e-01,2.985982741882e-02,1.685955262051e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,1,0,1),3, 1e-5);

    Sab  = {-3.861184580020e-02,-1.920105122015e-02,2.985982741882e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,1,1,0),3, 1e-5);

    Sab  = {-1.191723635809e-01,2.574123053346e-02,1.658879301045e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,0,2,0,0),3, 1e-5);

    Sab  = {-1.244159475784e-01,-1.343692233847e-01,-7.586798679230e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,0,0),3, 1e-5);

    Sab  = {-2.528932893077e-01,-2.731247524523e-01,-1.262473503268e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,0,1),3, 1e-5);

    Sab  = {-6.281844137199e-01,-6.784391668175e-01,-2.676045238273e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,0,2),3, 1e-5);

    Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,1,0),3, 1e-5);

    Sab  = {-9.104158415076e-02,4.074616677325e-01,-4.544904611764e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,1,1),3, 1e-5);

    Sab  = {-1.405402543846e-01,2.737548977757e-02,-8.570047788058e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,0,2,0),3, 1e-5);

    Sab  = {3.732478427352e-01,-4.478974112823e-02,-2.528932893077e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,1,0,0),3, 1e-5);

    Sab  = {7.586798679230e-01,-9.104158415075e-02,-4.208245010892e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,1,0,1),3, 1e-5);

    Sab  = {1.343692233847e-01,6.681965824611e-02,-9.104158415076e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,1,1,0),3, 1e-5);

    Sab  = {4.147198252614e-01,-8.957948225645e-02,-5.057865786153e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,0,1,2,0,0),3, 1e-5);

    Sab  = {-1.930592290010e-02,8.640473049066e-02,-1.343692233847e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,0,0),3, 1e-5);

    Sab  = {-4.478974112823e-02,2.004589747383e-01,-2.731247524523e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,0,1),3, 1e-5);

    Sab  = {-1.232181223176e-01,5.514695518836e-01,-6.784391668175e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,0,2),3, 1e-5);

    Sab  = {2.880157683022e-02,1.083293945770e-01,2.004589747383e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,1,0),3, 1e-5);

    Sab  = {6.681965824611e-02,2.513241954187e-01,4.074616677325e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,1,1),3, 1e-5);

    Sab  = {3.933260025514e-03,1.019700521751e-02,2.737548977757e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,0,2,0),3, 1e-5);

    Sab  = {5.791776870029e-02,2.880157683022e-02,-4.478974112823e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,1,0,0),3, 1e-5);

    Sab  = {1.343692233847e-01,6.681965824611e-02,-9.104158415076e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,1,0,1),3, 1e-5);

    Sab  = {-8.640473049066e-02,3.610979819234e-02,6.681965824611e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,1,1,0),3, 1e-5);

    Sab  = {6.435307633366e-02,5.760315366044e-02,-8.957948225645e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(1,1,0,2,0,0),3, 1e-5);

    Sab  = {-1.787585453713e-02,-8.365899923376e-02,-5.391357728398e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,0,0),3, 1e-5);

    Sab  = {-4.147198252614e-02,-1.940888782223e-01,-1.095870920333e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,0,1),3, 1e-5);

    Sab  = {-1.140908539978e-01,-5.339451967095e-01,-2.722132459453e+00};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,0,2),3, 1e-5);

    Sab  = {-6.435307633366e-03,1.248068329310e-01,-1.940888782223e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,1,0),3, 1e-5);

    Sab  = {-1.492991370941e-02,2.895518523998e-01,-3.945135313199e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,1,1),3, 1e-5);

    Sab  = {-2.019256528514e-02,1.704412677723e-02,-6.090077689998e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,0,2,0),3, 1e-5);

    Sab  = {2.085516362665e-01,-6.435307633366e-03,-4.147198252614e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,1,0,0),3, 1e-5);

    Sab  = {4.838397961383e-01,-1.492991370941e-02,-8.429776310255e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,1,0,1),3, 1e-5);

    Sab  = {7.507858905594e-02,9.600525610074e-03,-1.492991370941e-02};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,1,1,0),3, 1e-5);

    Sab  = {-1.191723635809e-02,-6.435307633366e-02,-4.147198252614e-01};
    CHECK_ARRAY_CLOSE(Sab, integrator.overlapIntegral_derivative(2,0,0,2,0,0),3, 1e-5);

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


    integrator.setMaxAngularMomentum(2);
    integrator.updateHermiteCoefficients(true, false);

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
TEST(GTOkineticIntegral_derivative)
{
    Integrator integrator;

    rowvec posA = {1.2,2.3,3.4};
    rowvec posB = {-1.3,1.4,-2.4};
    rowvec Tab  = {0, 0, 0};
    integrator.setCorePositionA(posA);
    integrator.setCorePositionB(posB);
    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);

    integrator.setMaxAngularMomentum(2);

    integrator.updateHermiteCoefficients(true, false);
    integrator.updateHermiteCoefficients_derivative(true,false);


    Tab  = {4.091139572785e-02,1.472810246203e-02,9.491443808862e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,0,0),3, 1e-5);

    Tab  = {5.510133486353e-02,1.983648055087e-02,1.005608330648e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,0,1),3, 1e-5);

    Tab  = {8.782064003279e-02,3.161543041181e-02,1.302754383914e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,0,2),3, 1e-5);

    Tab  = {8.550207133995e-03,-2.419618925033e-02,1.983648055087e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,1,0),3, 1e-5);

    Tab  = {5.503763389837e-03,-3.475286842201e-02,7.068592975089e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,1,1),3, 1e-5);

    Tab  = {5.320597488176e-02,7.753874778772e-03,1.234378617257e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,0,2,0),3, 1e-5);

    Tab  = {-3.523688446359e-03,8.550207133996e-03,5.510133486353e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,1,0,0),3, 1e-5);

    Tab  = {-2.144599160392e-02,5.503763389837e-03,1.963498048636e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,1,0,1),3, 1e-5);

    Tab  = {-3.327826283366e-03,-1.497968466466e-02,5.503763389837e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,1,1,0),3, 1e-5);

    Tab  = {2.727426381857e-02,2.121901115334e-02,1.367447385438e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,0,2,0,0),3, 1e-5);

    Tab  = {-8.265200229529e-02,-2.975472082630e-02,-1.508412495972e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,0,0),3, 1e-5);

    Tab  = {-2.945247072954e-02,-1.060288946263e-02,4.187293763460e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,0,1),3, 1e-5);

    Tab  = {1.061153609984e-01,3.820152995943e-02,3.732782385369e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,0,2),3, 1e-5);

    Tab  = {-8.255645084755e-03,5.212930263301e-02,-1.060288946264e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,1,0),3, 1e-5);

    Tab  = {3.309707092973e-02,3.154992602106e-02,8.779273133667e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,1,1),3, 1e-5);

    Tab  = {-1.176974704839e-01,-3.136356259454e-02,-2.198521566410e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,0,2,0),3, 1e-5);

    Tab  = {3.216898740587e-02,-8.255645084755e-03,-2.945247072955e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,1,0,0),3, 1e-5);

    Tab  = {1.115712886245e-01,3.309707092973e-02,2.438686981574e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,1,0,1),3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092972e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,1,1,0),3, 1e-5);

    Tab  = {-5.510133486352e-02,-3.084400733054e-02,-1.398307940375e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,1,2,0,0),3, 1e-5);

    Tab  = {1.464571954140e-01,5.272459034902e-02,1.744766887699e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,0,0),3, 1e-5);

    Tab  = {-2.280497100773e-01,-8.209789562783e-02,-6.856183990998e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,0,1),3, 1e-5);

    Tab  = {-9.826166149089e-01,-3.537419813672e-01,-1.763373539942e+00};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,0,2),3, 1e-5);

    Tab  = {-2.827008855960e-02,-1.078153621574e-01,-8.209789562774e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,1,0),3, 1e-5);

    Tab  = {-2.270073992127e-01,7.031047633496e-02,-4.416162986024e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,1,1),3, 1e-5);

    Tab  = {2.571120252525e-01,1.302537805037e-01,3.611049576180e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,0,2,0),3, 1e-5);

    Tab  = {-1.761661540526e-01,-2.827008855960e-02,-2.280497100771e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,1,0,0),3, 1e-5);

    Tab  = {-4.785429688727e-01,-2.270073992127e-01,-1.226711940562e+00};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,1,0,1),3, 1e-5);

    Tab  = {-9.041804176183e-02,1.301669956255e-02,-2.270073992127e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,1,1,0),3, 1e-5);

    Tab  = {9.763813027596e-02,-2.543724513444e-03,-1.877484875895e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,0,2,2,0,0),3, 1e-5);

    Tab  = {-1.282531070099e-02,3.629428387549e-02,-2.975472082630e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,0,0),3, 1e-5);

    Tab  = {-8.255645084755e-03,5.212930263301e-02,-1.060288946263e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,0,1),3, 1e-5);

    Tab  = {1.172155309473e-02,9.204039914690e-02,3.820152995943e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,0,2),3, 1e-5);

    Tab  = {2.246952699699e-02,2.518944398691e-02,5.212930263301e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,1,0),3, 1e-5);

    Tab  = {2.005586667488e-02,1.822763878263e-02,3.154992602106e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,1,1),3, 1e-5);

    Tab  = {-1.351877698041e-02,1.837984583950e-02,-3.136356259454e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,0,2,0),3, 1e-5);

    Tab  = {4.991739425050e-03,2.246952699699e-02,-8.255645084754e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,1,0,0),3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092973e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,1,0,1),3, 1e-5);

    Tab  = {-6.334914546175e-03,7.856740854583e-03,2.005586667488e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,1,1,0),3, 1e-5);

    Tab  = {-8.550207133995e-03,5.415555857965e-02,-3.084400733054e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,0,2,0,0),3, 1e-5);

    Tab  = {1.238346762713e-02,-7.819395394952e-02,1.590433419396e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,0,0),3, 1e-5);

    Tab  = {-4.964560639459e-02,-4.732488903159e-02,-1.316890970050e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,0,1),3, 1e-5);

    Tab  = {-2.201276949753e-01,2.686939080730e-02,-4.327805573967e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,0,2),3, 1e-5);

    Tab  = {-3.008380001232e-02,-2.734145817395e-02,-4.732488903160e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,1,0),3, 1e-5);

    Tab  = {5.046591122420e-02,8.436186990016e-02,1.571926473900e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,1,1),3, 1e-5);

    Tab  = {5.361546143843e-02,-5.828417101632e-02,1.108690935568e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,0,2,0),3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,1,0,0),3, 1e-5);

    Tab  = {-8.209847605341e-02,5.046591122421e-02,-2.407668076873e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,1,0,1),3, 1e-5);

    Tab  = {3.808222059976e-02,3.297634872656e-02,5.046591122419e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,1,1,0),3, 1e-5);

    Tab  = {8.255645084751e-03,-9.224103598278e-02,-5.559125239680e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,1,1,2,0,0),3, 1e-5);

    Tab  = {6.857419882415e-02,-9.639098252967e-04,1.590921412720e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,0,0),3, 1e-5);

    Tab  = {1.076695371465e-01,2.224974320323e-02,2.040771936304e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,0,1),3, 1e-5);

    Tab  = {2.138436938137e-01,1.004268359624e-01,3.525579867858e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,0,2),3, 1e-5);

    Tab  = {9.590406553117e-03,2.675467803668e-03,2.224974320323e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,1,0),3, 1e-5);

    Tab  = {4.673575720113e-02,-1.484308548881e-02,1.020333523379e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,1,1),3, 1e-5);

    Tab  = {1.078924496974e-01,-9.834808072682e-04,2.503104832979e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,0,2,0),3, 1e-5);

    Tab  = {6.931507033685e-04,9.590406553116e-03,1.076695371465e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,1,0,0),3, 1e-5);

    Tab  = {-1.553275841004e-02,4.673575720113e-02,9.955336244010e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,1,0,1),3, 1e-5);

    Tab  = {1.375111856278e-02,-6.397881676221e-03,4.673575720115e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,1,1,0),3, 1e-5);

    Tab  = {4.571613254943e-02,1.214460218731e-02,2.496208103767e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(0,2,0,2,0,0),3, 1e-5);

    Tab  = {5.285532669538e-03,-1.282531070099e-02,-8.265200229529e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,0,0),3, 1e-5);

    Tab  = {3.216898740587e-02,-8.255645084755e-03,-2.945247072954e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,0,1),3, 1e-5);

    Tab  = {1.203805097404e-01,1.172155309473e-02,1.061153609984e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,0,2),3, 1e-5);

    Tab  = {4.991739425050e-03,2.246952699699e-02,-8.255645084755e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,1,0),3, 1e-5);

    Tab  = {1.874719404663e-02,2.005586667488e-02,3.309707092973e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,1,1),3, 1e-5);

    Tab  = {2.474306569721e-03,-1.351877698041e-02,-1.176974704839e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,0,2,0),3, 1e-5);

    Tab  = {6.136709359178e-02,4.991739425050e-03,3.216898740587e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,1,0,0),3, 1e-5);

    Tab  = {8.265200229529e-02,1.874719404663e-02,1.115712886245e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,1,0,1),3, 1e-5);

    Tab  = {1.282531070099e-02,-6.334914546176e-03,1.874719404663e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,1,1,0),3, 1e-5);

    Tab  = {1.670319847949e-02,-8.550207133995e-03,-5.510133486352e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,0,2,0,0),3, 1e-5);

    Tab  = {-4.825348110881e-02,1.238346762713e-02,4.417870609432e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,0,0),3, 1e-5);

    Tab  = {-1.673569329367e-01,-4.964560639459e-02,-3.658030472361e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,0,1),3, 1e-5);

    Tab  = {-5.053504583775e-01,-2.201276949754e-01,-1.202168214991e+00};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,0,2),3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,1,0),3, 1e-5);

    Tab  = {-8.209847605341e-02,5.046591122420e-02,-2.407668076873e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,1,1),3, 1e-5);

    Tab  = {-4.234024791492e-02,5.361546143843e-02,1.240970880481e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,0,2,0),3, 1e-5);

    Tab  = {-1.239780034429e-01,-2.812079106995e-02,-1.673569329367e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,1,0,0),3, 1e-5);

    Tab  = {-4.417870609430e-02,-8.209847605341e-02,-4.249279898629e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,1,0,1),3, 1e-5);

    Tab  = {-1.238346762713e-02,3.808222059976e-02,-8.209847605341e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,1,1,0),3, 1e-5);

    Tab  = {4.140562735410e-02,8.255645084759e-03,2.945247072955e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,0,1,2,0,0),3, 1e-5);

    Tab  = {-7.487609137574e-03,-3.370429049548e-02,1.238346762713e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,0,0),3, 1e-5);

    Tab  = {-2.812079106995e-02,-3.008380001232e-02,-4.964560639459e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,0,1),3, 1e-5);

    Tab  = {-9.457782461603e-02,-5.707906268285e-03,-2.201276949753e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,0,2),3, 1e-5);

    Tab  = {9.502371819266e-03,-1.178511128187e-02,-3.008380001232e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,1,0),3, 1e-5);

    Tab  = {3.808222059976e-02,3.297634872655e-02,5.046591122421e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,1,1),3, 1e-5);

    Tab  = {9.591335708574e-03,-2.512248750705e-02,5.361546143843e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,0,2,0),3, 1e-5);

    Tab  = {-1.923796605149e-02,9.502371819265e-03,-2.812079106995e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,1,0,0),3, 1e-5);

    Tab  = {-1.238346762712e-02,3.808222059975e-02,-8.209847605341e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,1,0,1),3, 1e-5);

    Tab  = {3.370429049548e-02,2.207068427119e-02,3.808222059976e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,1,1,0),3, 1e-5);

    Tab  = {6.425011141151e-03,-2.246952699699e-02,8.255645084759e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(1,1,0,2,0,0),3, 1e-5);

    Tab  = {1.022784893196e-02,2.933264701749e-02,1.890326141127e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,0,0),3, 1e-5);

    Tab  = {1.377533371588e-02,2.147041030723e-02,8.404514972528e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,0,1),3, 1e-5);

    Tab  = {2.195516000820e-02,-1.553924858651e-02,-1.796618623990e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,0,2),3, 1e-5);

    Tab  = {2.137551783499e-03,-5.098810130656e-02,2.147041030723e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,1,0),3, 1e-5);

    Tab  = {1.375940847460e-03,-4.879995045527e-02,-6.442699361569e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,1,1),3, 1e-5);

    Tab  = {1.330149372044e-02,2.897602265551e-02,2.662544063993e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,0,2,0),3, 1e-5);

    Tab  = {-2.065018716129e-02,2.137551783499e-03,1.377533371588e-02};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,1,0,0),3, 1e-5);

    Tab  = {6.821311685898e-02,1.375940847467e-03,4.908745121587e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,1,0,1),3, 1e-5);

    Tab  = {1.058479399536e-02,-3.744921166164e-03,1.375940847460e-03};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,1,1,0),3, 1e-5);

    Tab  = {1.097834880885e-01,5.947253902451e-02,3.832674737135e-01};
    CHECK_ARRAY_CLOSE(Tab, integrator.kineticIntegral_derivative(2,0,0,2,0,0),3, 1e-5);




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


    integrator.setExponentA(0.2);
    integrator.setExponentB(0.3);
    integrator.setMaxAngularMomentum(2);
    integrator.updateHermiteCoefficients(true, false);


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
