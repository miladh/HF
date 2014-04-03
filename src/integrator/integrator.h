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

    rowvec corePositionA() const;
    void setCorePositionA(const rowvec &corePositionA);

    rowvec corePositionB() const;
    void setCorePositionB(const rowvec &corePositionB);

    rowvec corePositionC() const;
    void setCorePositionC(const rowvec &corePositionC);

    rowvec corePositionD() const;
    void setCorePositionD(const rowvec &corePositionD);


    double kineticIntegral();
    double nuclearAttractionIntegral();
    double electronRepulsionIntegral();

    void updateHermiteCoefficients(bool oneParticleIntegral, bool twoParticleIntegral, bool kin= true);
    void updateHermiteCoefficients_derivative(bool oneParticleIntegral, bool twoParticleIntegral, bool kin =true);


   rowvec overlapIntegral_derivative();
   rowvec kineticIntegral_derivative();
   rowvec nuclearAttractionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                               bool differentiateWrtC);

   rowvec electronRepulsionIntegral_derivative(bool differentiateWrtA, bool differentiateWrtB,
                                               bool differentiateWrtC,  bool differentiateWrtD);

   double overlapIntegral();


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





    double overlapIntegral(int cor, int iA, int iB);
    double kineticIntegral(int cor, int iA, int iB);
    double overlapIntegral_derivative(int cor, int iA, int iB);
    double kineticIntegral_derivative(int cor, int iA, int iB);

    rowvec nuclearAttractionIntegral_R_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_P_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec nuclearAttractionIntegral_C_derivative(int iA, int jA, int kA, int iB, int jB, int kB);
    rowvec electronRepulsionIntegral_Rab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);
    rowvec electronRepulsionIntegral_Rcd_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);
    rowvec electronRepulsionIntegral_Pab_derivative(int iA, int jA, int kA, int iB, int jB, int kB,
                                                  int iC, int jC, int kC, int iD, int jD, int kD);


    rowvec electronRepulsionIntegral_Pcd_derivative(int iA, int jA, int kA, int iB, int jB, int kB, int iC, int jC, int kC, int iD, int jD, int kD);
};
}
#endif // INTEGRATOR_H

