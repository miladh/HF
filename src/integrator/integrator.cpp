#include "integrator.h"

using namespace hf;

Integrator::Integrator(const int maxAngularMomentum):
    m_nuclearSourceCharge(rowvec(3))
{
    Eab = new HermiteCoefficients(maxAngularMomentum);
    Ecd = new HermiteCoefficients(maxAngularMomentum);

    m_overlap = new OverlapIntegral(Eab->coefficients(), &m_primitiveA, &m_primitiveB);
    m_dipole = new DipoleIntegral(Eab->coefficients(), &m_primitiveA, &m_primitiveB);
    m_kinetic = new KineticIntegral(m_overlap, &m_primitiveA, &m_primitiveB);
    m_nuclearAttraction = new NuclearAttractionIntegral(2 * maxAngularMomentum + 1,
                                                        Eab->coefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_nuclearSourceCharge);

    m_electronRepulsion = new ElectronRepulsionIntegral(4 * maxAngularMomentum + 1,
                                                        Eab->coefficients(),
                                                        Ecd->coefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_primitiveC,
                                                        &m_primitiveD);



    m_overlapGD = new OverlapIntegralGD(m_overlap, Eab->QDerivativeCoefficients());
    m_kineticGD = new KineticIntegralGD(m_kinetic, m_overlapGD);
    m_nuclearAttractionGD = new NuclearAttractionIntegralGD(2 * maxAngularMomentum + 2,
                                                        Eab->coefficients(),
                                                        Eab->QDerivativeCoefficients(),
                                                        &m_primitiveA,
                                                        &m_primitiveB,
                                                        &m_nuclearSourceCharge);

    m_electronRepulsionGD = new ElectronRepulsionIntegralGD(4 * maxAngularMomentum + 3,
                                                            Eab->coefficients(),
                                                            Ecd->coefficients(),
                                                            Eab->QDerivativeCoefficients(),
                                                            Ecd->QDerivativeCoefficients(),
                                                            &m_primitiveA,
                                                            &m_primitiveB,
                                                            &m_primitiveC,
                                                            &m_primitiveD);

}

void Integrator::setNuclearSourceCharge(const rowvec &nuclearSourceCharge)
{
    m_nuclearSourceCharge = nuclearSourceCharge;
}

void Integrator::setPrimitiveA(const PrimitiveGTO &primitiveA)
{
    m_primitiveA = primitiveA;
}

void Integrator::setPrimitiveB(const PrimitiveGTO &primitiveB)
{
    m_primitiveB = primitiveB;
}

void Integrator::setPrimitiveC(const PrimitiveGTO &primitiveC)
{
    m_primitiveC = primitiveC;
}

void Integrator::setPrimitiveD(const PrimitiveGTO &primitiveD)
{
    m_primitiveD = primitiveD;
}

/********************************************************************************************
 *
 *                                Update functions
 *
 * ******************************************************************************************/
void Integrator::updateKineticHermiteCoefficients()
{
    Eab->updateE(m_primitiveA, m_primitiveB);
}

void Integrator::updateOverlapHermiteCoefficients()
{
    Eab->updateE(m_primitiveA, m_primitiveB,false);
}

void Integrator::updateElectronRepulsionHermiteCoefficients()
{
    Ecd->updateE(m_primitiveC, m_primitiveD, false);
}

void Integrator::updateNuclearAttractionHermiteIntegrals()
{
    m_nuclearAttractionGD->updateHermiteIntegrals();
}

void Integrator::updateElectronRepulsionHermiteIntegrals()
{
    m_electronRepulsionGD->updateHermiteIntegrals();
}


void Integrator::updateKineticHermiteCoefficientsGD()
{
    Eab->updatedE_dQ(m_primitiveA, m_primitiveB);
}

void Integrator::updateOverlapHermiteCoefficientsGD()
{
    Eab->updatedE_dQ(m_primitiveA, m_primitiveB,false);
}

void Integrator::updateElectronRepulsionHermiteCoefficientsGD()
{
    Ecd->updatedE_dQ(m_primitiveC, m_primitiveD, false);
}
/********************************************************************************************
 *
 *                              Molecular Gaussian Integrals
 *
 * ******************************************************************************************/
double Integrator::overlapIntegral()
{
    return m_overlap->evaluate();
}

rowvec Integrator::dipoleIntegral()
{
    return m_dipole->evaluate();
}

double Integrator::kineticIntegral()
{
    return m_kinetic->evaluate();
}

double Integrator::nuclearAttractionIntegral()
{
    return m_nuclearAttraction->evaluate();
}

double Integrator::electronRepulsionIntegral()
{
    return m_electronRepulsion->evaluate();
}


/********************************************************************************************
 *
 *                  Molecular Gaussian Integral Geometrical Derivatives (GD)
 *
 * ******************************************************************************************/

rowvec Integrator::QDerivativeOverlapIntegral()
{
    return m_overlapGD->evaluate();
}


rowvec Integrator::QDerivativeKineticIntegral() {

    return m_kineticGD->evaluate();
}


rowvec Integrator::QDerivativeNuclearAttractionIntegral()
{

    return m_nuclearAttractionGD->QDerivative();
}

rowvec Integrator::PDerivativeNuclearAttractionIntegral()
{

    return m_nuclearAttractionGD->PDerivative();
}


rowvec Integrator::QabDerivativeElectronRepulsionIntegral()
{

    return m_electronRepulsionGD->QabDerivative();
}

rowvec Integrator::PabDerivativeElectronRepulsionIntegral()
{

    return m_electronRepulsionGD->PabDerivative();
}

rowvec Integrator::QcdDerivativeElectronRepulsionIntegral()
{

    return m_electronRepulsionGD->QcdDerivative();
}

rowvec Integrator::PcdDerivativeElectronRepulsionIntegral()
{

    return m_electronRepulsionGD->PcdDerivative();
}












