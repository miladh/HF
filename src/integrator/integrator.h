#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <iostream>
#include <armadillo>

#include "../math/hermitecoefficients.h"
#include "../math/hermiteintegrals.h"
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
    Integrator();

    int maxAngularMomentum() const;
    void setMaxAngularMomentum(const int maxAngularMomentum);


    void setPrimitiveA(const PrimitiveGTO &primitiveA);
    void setPrimitiveB(const PrimitiveGTO &primitiveB);
    void setPrimitiveC(const PrimitiveGTO &primitiveC);
    void setPrimitiveD(const PrimitiveGTO &primitiveD);

    rowvec corePositionC() const;
    void setCorePositionC(const rowvec &corePositionC);

    double overlapIntegral();
    double kineticIntegral();
    double nuclearAttractionIntegral();
    double electronRepulsionIntegral();


    void updateHermiteIntegrals();
    void updateHermiteCoefficients(bool oneParticleIntegral, bool twoParticleIntegral, bool kin= true);
    void updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral, bool kin =true);


   rowvec QDerivativeOverlapIntegral();
   rowvec QDerivativeKineticIntegral();

   rowvec QDerivativeNuclearAttractionIntegral();
   rowvec PDerivativeNuclearAttractionIntegral();
   rowvec CDerivativeNuclearAttractionIntegral();

   rowvec QabDerivativeElectronRepulsionIntegral();
   rowvec PabDerivativeElectronRepulsionIntegral();
   rowvec QcdDerivativeElectronRepulsionIntegral();
   rowvec PcdDerivativeElectronRepulsionIntegral();




   rowvec nuclearAttractionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                               bool differentiateWrtC);

   rowvec electronRepulsionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                               bool differentiateWrtC,  bool differentiateWrtD);


private:
   PrimitiveGTO m_primitiveA;
   PrimitiveGTO m_primitiveB;
   PrimitiveGTO m_primitiveC;
   PrimitiveGTO m_primitiveD;

   rowvec m_corePositionC;

   field<cube> m_Eab, m_Ecd;
   field<cube> m_dEab, m_dEcd;
   field<cube> m_Ree, m_Ren;

    HermiteCoefficients m_hermiteCoefficients;
    HermiteIntegrals *m_hermiteIntegrals;

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



    rowvec nuclearAttractionIntegral_R_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_P_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_C_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec electronRepulsionIntegral_Rab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);
    rowvec electronRepulsionIntegral_Rcd_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);
    rowvec electronRepulsionIntegral_Pab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);


    rowvec electronRepulsionIntegral_Pcd_derivative(int iA, int jA, int kA,
                                                    int iB, int jB, int kB,
                                                    int iC, int jC, int kC,
                                                    int iD, int jD, int kD);
};
}
#endif // INTEGRATOR_H

