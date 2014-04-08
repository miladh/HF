#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>

#include "../math/hermitecoefficients.h"
#include "../primitiveGTO/primitiveGTO.h"

#include "overlap/overlapintegral.h"
#include "kinetic/kineticintegral.h"
#include "nuclearAttraction/nuclearattractionintegral.h"
#include "electronRepulsion/electronrepulsionintegral.h"

#include "overlap/overlapintegralgd.h"
#include "kinetic/kineticintegralgd.h"
#include "nuclearAttraction/nuclearattractionintegralgd.h"
#include "electronRepulsion/electronrepulsionintegralgd.h"



using namespace arma;
using namespace std;

namespace hf
{
class Integrator
{
public:
    Integrator(const int maxAngularMomentum);

    void setPrimitiveA(const PrimitiveGTO &primitiveA);
    void setPrimitiveB(const PrimitiveGTO &primitiveB);
    void setPrimitiveC(const PrimitiveGTO &primitiveC);
    void setPrimitiveD(const PrimitiveGTO &primitiveD);

    rowvec nuclearSourceCharge() const;
    void setNuclearSourceCharge(const rowvec &nuclearSourceCharge);

    void updateNuclearAttractionHermiteIntegrals();
    void updateElectronRepulsionHermiteIntegrals();


    double overlapIntegral();
    double kineticIntegral();
    double nuclearAttractionIntegral();
    double electronRepulsionIntegral();

   rowvec QDerivativeOverlapIntegral();
   rowvec QDerivativeKineticIntegral();

   rowvec QDerivativeNuclearAttractionIntegral();
   rowvec PDerivativeNuclearAttractionIntegral();
   rowvec CDerivativeNuclearAttractionIntegral();

   rowvec QabDerivativeElectronRepulsionIntegral();
   rowvec PabDerivativeElectronRepulsionIntegral();
   rowvec QcdDerivativeElectronRepulsionIntegral();
   rowvec PcdDerivativeElectronRepulsionIntegral();

   void updateOverlapHermiteCoefficients();
   void updateKineticHermiteCoefficients();
   void updateElectronRepulsionHermiteCoefficientsGD();

   void updateOverlapHermiteCoefficientsGD();
   void updateKineticHermiteCoefficientsGD();
   void updateElectronRepulsionHermiteCoefficients();


private:
   PrimitiveGTO m_primitiveA;
   PrimitiveGTO m_primitiveB;
   PrimitiveGTO m_primitiveC;
   PrimitiveGTO m_primitiveD;

   rowvec m_nuclearSourceCharge;

    HermiteCoefficients* Eab;
    HermiteCoefficients* Ecd;

    OverlapIntegral* m_overlap;
    KineticIntegral* m_kinetic;
    NuclearAttractionIntegral* m_nuclearAttraction;
    ElectronRepulsionIntegral* m_electronRepulsion;

    OverlapIntegralGD* m_overlapGD;
    KineticIntegralGD* m_kineticGD;
    NuclearAttractionIntegralGD* m_nuclearAttractionGD;
    ElectronRepulsionIntegralGD* m_electronRepulsionGD;

};
}
#endif // INTEGRATOR_H

