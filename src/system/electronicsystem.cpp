#include "electronicsystem.h"

using namespace hf;

ElectronicSystem::ElectronicSystem()
{
}

void ElectronicSystem::addAtoms(vector<Atom*> atoms)
{
    int maxAngularMomentum = 0;
    int id = 0;
    m_atoms = atoms;
    for(Atom* atom : m_atoms){
        m_nElectrons        += atom->nElectrons();
        maxAngularMomentum = max(maxAngularMomentum , atom->angularMomentum());

        for(const ContractedGTO &CGTO : atom->contractedGTOs()){
            m_basisFunctions.push_back(&CGTO);
            m_basisFucntionIndexToAtomID.push_back(id);
        }
        id++;
    }


    m_nAtoms             = m_atoms.size();
    m_nBasisFunctions    = m_basisFunctions.size();
    m_nSpinDownElectrons = ceil(m_nElectrons/2.0);
    m_nSpinUpElectrons   = floor(m_nElectrons/2.0);

    if(m_nSpinDownElectrons + m_nSpinUpElectrons !=m_nElectrons){
        throw logic_error("Number of electrons not conserved in system!");
    }

    integrator = new Integrator(maxAngularMomentum);
}

/********************************************************************************************
 *
 *                                  System properties functions
 *
 * ******************************************************************************************/
const int& ElectronicSystem::nElectrons() const
{
    return m_nElectrons;
}

const int& ElectronicSystem::nSpinUpElectrons() const
{
    return m_nSpinUpElectrons;
}

const int& ElectronicSystem::nSpinDownElectrons() const
{
    return m_nSpinDownElectrons;
}

const int& ElectronicSystem::nAtoms()
{
    return m_nAtoms;
}

const int& ElectronicSystem::nBasisFunctions()
{
    return m_nBasisFunctions;
}

vector<const ContractedGTO *> ElectronicSystem::basisFunctions() const
{
    return m_basisFunctions;
}
vector<Atom *> ElectronicSystem::atoms() const
{
    return m_atoms;
}

/********************************************************************************************
 *
 *                                  Molecular Integrals
 *
 * ******************************************************************************************/

double ElectronicSystem::overlapIntegral(const int& p, const int& q)
{
    double Spq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateOverlapHermiteCoefficients();

            Spq += integrator->overlapIntegral();
        }
    }
    return Spq;
}

double ElectronicSystem::oneParticleIntegral(const int& p, const int& q)
{
    double hpq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateKineticHermiteCoefficients();

            hpq += integrator->kineticIntegral();

            for(const Atom *atom : m_atoms){
                integrator->setNuclearSourceCharge(atom->corePosition());
                hpq -= atom->coreCharge() * integrator->nuclearAttractionIntegral();

            }

        }
    }

    return hpq;
}


double ElectronicSystem::twoParticleIntegral(const int& p, const int& q,
                                             const int& r, const int& s)
{
    double Qpqrs = 0.0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    const ContractedGTO *CGr = m_basisFunctions.at(r);
    const ContractedGTO *CGs = m_basisFunctions.at(s);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateOverlapHermiteCoefficients();


            for(const PrimitiveGTO &Gc : CGr->primitiveGTOs()) {
                integrator->setPrimitiveC(Gc);

                for(const PrimitiveGTO &Gd : CGs->primitiveGTOs()) {
                    integrator->setPrimitiveD(Gd);
                    integrator->updateElectronRepulsionHermiteCoefficients();

                    Qpqrs += integrator->electronRepulsionIntegral();

                }
            }
        }
    }

    return Qpqrs;
}


double ElectronicSystem::nuclearPotential()
{
    double Vn = 0;

    for(int i = 0; i < int(m_atoms.size()); i++){
        for(int j = i+1; j < int(m_atoms.size()); j++){
            rowvec AB = m_atoms.at(i)->corePosition() - m_atoms.at(j)->corePosition();

            Vn += m_atoms.at(i)->coreCharge() * m_atoms.at(j)->coreCharge() / sqrt(dot(AB,AB));
        }
    }

    return Vn;
}

/********************************************************************************************
 *
 *                        Molecular Integral Geometrical Derivatives (GD)
 *
 * ******************************************************************************************/

mat ElectronicSystem::nuclearPotentialGD()
{
    mat dVn = zeros(m_nAtoms, 3);
    for(int i = 0; i < m_nAtoms; i++){
        for(int j = i+1; j < m_nAtoms; j++){
            rowvec AB = m_atoms.at(i)->corePosition() - m_atoms.at(j)->corePosition();
            rowvec dV = AB * m_atoms.at(i)->coreCharge() * m_atoms.at(j)->coreCharge()/pow(dot(AB,AB),1.5);

            dVn.row(i) -= dV;
            dVn.row(j) += dV;

        }
    }

    return dVn;
}

mat ElectronicSystem::overlapIntegralGD(const int& q, const int& p)
{
    mat dSpq = zeros(m_nAtoms, 3);

    int coreA = m_basisFucntionIndexToAtomID.at(p);
    int coreB = m_basisFucntionIndexToAtomID.at(q);

    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateOverlapHermiteCoefficients();
            integrator->updateOverlapHermiteCoefficientsGD();



            dSpq.row(coreA) += integrator->QDerivativeOverlapIntegral();

        }
    }

    dSpq.row(coreB) -= dSpq.row(coreA);

    return dSpq;
}

mat ElectronicSystem::oneParticleIntegralGD(const int& q, const int& p)
{
    mat dhpq = zeros(m_nAtoms, 3);

    int coreA = m_basisFucntionIndexToAtomID.at(p);
    int coreB = m_basisFucntionIndexToAtomID.at(q);

    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateKineticHermiteCoefficients();
            integrator->updateKineticHermiteCoefficientsGD();

            rowvec QDerivative = integrator->QDerivativeKineticIntegral();
            dhpq.row(coreA) += QDerivative;
            dhpq.row(coreB) -= QDerivative;

            double ab   = Ga.exponent() + Gb.exponent();
            double a_ab = Ga.exponent()/ab;
            double b_ab = Gb.exponent()/ab;
            for(int i = 0; i < m_nAtoms; i++){
                integrator->setNuclearSourceCharge(m_atoms.at(i)->corePosition());
                integrator->updateNuclearAttractionHermiteIntegrals();

                rowvec QDerivative = integrator->QDerivativeNuclearAttractionIntegral();
                rowvec PDerivative = integrator->PDerivativeNuclearAttractionIntegral();

                dhpq.row(coreA) -= (a_ab * PDerivative + QDerivative) * m_atoms.at(i)->coreCharge() ;
                dhpq.row(coreB) -= (b_ab * PDerivative - QDerivative) * m_atoms.at(i)->coreCharge() ;

                dhpq.row(i) +=m_atoms.at(i)->coreCharge()
                        * integrator->CDerivativeNuclearAttractionIntegral();

            }

        }
    }

    return dhpq;
}


mat ElectronicSystem::twoParticleIntegralGD(const int& p, const int& q,
                                             const int& r, const int& s)
{
    mat dQpqrs = zeros(m_nAtoms, 3);

    int coreA = m_basisFucntionIndexToAtomID.at(p);
    int coreB = m_basisFucntionIndexToAtomID.at(q);
    int coreC = m_basisFucntionIndexToAtomID.at(r);
    int coreD = m_basisFucntionIndexToAtomID.at(s);

    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    const ContractedGTO *CGr = m_basisFunctions.at(r);
    const ContractedGTO *CGs = m_basisFunctions.at(s);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            integrator->setPrimitiveB(Gb);
            integrator->updateOverlapHermiteCoefficients();
            integrator->updateOverlapHermiteCoefficientsGD();


            double ab   = Ga.exponent() + Gb.exponent();
            double a_ab = Ga.exponent()/ab;
            double b_ab = Gb.exponent()/ab;
            for(const PrimitiveGTO &Gc : CGr->primitiveGTOs()) {
                integrator->setPrimitiveC(Gc);

                for(const PrimitiveGTO &Gd : CGs->primitiveGTOs()) {
                    integrator->setPrimitiveD(Gd);
                    integrator->updateElectronRepulsionHermiteCoefficients();
                    integrator->updateElectronRepulsionHermiteCoefficientsGD();


                    integrator->updateElectronRepulsionHermiteIntegrals();
                    double cd   = Gc.exponent() + Gd.exponent();
                    double c_cd = Gc.exponent()/cd;
                    double d_cd = Gd.exponent()/cd;

                    double normalizationFactor = 2.0 *pow(M_PI,2.5)/ (ab*cd*sqrt(ab+cd));

                    rowvec QabDerivative = integrator->QabDerivativeElectronRepulsionIntegral();
                    rowvec PabDerivative = integrator->PabDerivativeElectronRepulsionIntegral();
                    rowvec QcdDerivative = integrator->QcdDerivativeElectronRepulsionIntegral();
                    rowvec PcdDerivative = integrator->PcdDerivativeElectronRepulsionIntegral();

                    dQpqrs.row(coreA) += (a_ab * PabDerivative + QabDerivative) * normalizationFactor;
                    dQpqrs.row(coreB) += (b_ab * PabDerivative - QabDerivative) * normalizationFactor;
                    dQpqrs.row(coreC) += (c_cd * PcdDerivative + QcdDerivative) * normalizationFactor;
                    dQpqrs.row(coreD) += (d_cd * PcdDerivative - QcdDerivative) * normalizationFactor;


                }
            }
        }
    }

    return dQpqrs;
}



















