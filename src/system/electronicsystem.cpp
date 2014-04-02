#include "electronicsystem.h"

using namespace hf;

ElectronicSystem::ElectronicSystem()
{
}

void ElectronicSystem::addAtoms(vector<Atom*> atoms)
{
    int maxAngularMomentum = 0;

    m_atoms = atoms;
    for(Atom* atom : m_atoms){
        m_nElectrons        += atom->nElectrons();
        maxAngularMomentum = max(maxAngularMomentum , atom->angularMomentum());

        for(const ContractedGTO &CGTO : atom->contractedGTOs()){
            m_basisFunctions.push_back(&CGTO);
        }
    }

    integrator.setMaxAngularMomentum(maxAngularMomentum);

    m_nAtoms             = m_atoms.size();
    m_nBasisFunctions    = m_basisFunctions.size();
    m_nSpinDownElectrons = ceil(m_nElectrons/2.0);
    m_nSpinUpElectrons   = floor(m_nElectrons/2.0);

    if(m_nSpinDownElectrons + m_nSpinUpElectrons !=m_nElectrons){
        throw logic_error("Number of electrons not conserved in system!");
    }
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

    for(const PrimitiveGTO &Gp : CGp->primitiveGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitiveGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false, false);

            Spq += integrator.overlapIntegral();
        }
    }
    return Spq;
}

double ElectronicSystem::oneParticleIntegral(const int& p, const int& q)
{
    double hpq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);

    for(const PrimitiveGTO &Gp : CGp->primitiveGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitiveGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false,true);

            hpq += integrator.kineticIntegral();

            for(const Atom *atom : m_atoms){
                integrator.setCorePositionC(atom->corePosition());
                hpq -= atom->coreCharge() * integrator.nuclearAttractionIntegral();

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

    for(const PrimitiveGTO &Gp : CGp->primitiveGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitiveGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false);

            for(const PrimitiveGTO &Gr : CGr->primitiveGTOs()) {
                integrator.setPrimitiveC(Gr);

                for(const PrimitiveGTO &Gs : CGs->primitiveGTOs()) {
                    integrator.setPrimitiveD(Gs);
                    integrator.updateHermiteCoefficients(false, true);

                    Qpqrs += integrator.electronRepulsionIntegral();

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

rowvec ElectronicSystem::nuclearPotentialGD(int activeCore)
{

    rowvec dVn = {0,0,0};

    for(int i = 0; i < int(m_atoms.size()); i++){
        for(int j = i+1; j < int(m_atoms.size()); j++){
            rowvec AB = m_atoms.at(i)->corePosition() - m_atoms.at(j)->corePosition();

            if(activeCore == i){
                dVn -= AB * m_atoms.at(i)->coreCharge() * m_atoms.at(j)->coreCharge()/pow(dot(AB,AB),1.5);
            }else if(activeCore == j){
                dVn += AB * m_atoms.at(i)->coreCharge() * m_atoms.at(j)->coreCharge()/pow(dot(AB,AB),1.5);
            }
        }
    }

    return dVn;
}


mat ElectronicSystem::getOneParticleDerivative(const int a, const int b, const int N)
{
    mat dhab = zeros(2,3);
//    bool differentiateWrtA;
//    bool differentiateWrtB;
//    bool differentiateWrtC;

//    if(N == m_coreID.at(a)){
//        differentiateWrtA = true;
//    }else{
//        differentiateWrtA = false;
//    }

//    if(N == m_coreID.at(b)){
//        differentiateWrtB = true;
//    }else{
//        differentiateWrtB = false;
//    }


    //    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    //    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    //    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    //    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));

    //    integrator.setCorePositionA(coreA->corePosition());
    //    integrator.setCorePositionB(coreB->corePosition());

    //    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
    //        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
    //        integrator.setPrimitiveA(primitiveA);

    //        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
    //            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
    //            integrator.setPrimitiveB(primitiveB);

    //            integrator.updateHermiteCoefficients(true, false);
    //            integrator.updateHermiteCoefficients_derivative(true,false);


    //            if(differentiateWrtA){
    //                dhab.row(0) += integrator.overlapIntegral_derivative()
    //                        * primitiveA.weight() * primitiveB.weight();

    //                dhab.row(1) += primitiveA.weight() * primitiveB.weight() *
    //                        integrator.kineticIntegral_derivative();
    //            }
    //            else if(differentiateWrtB){
    //                dhab.row(0) -= integrator.overlapIntegral_derivative()
    //                        * primitiveA.weight() * primitiveB.weight();

    //                dhab.row(1) -= primitiveA.weight() * primitiveB.weight() *
    //                        integrator.kineticIntegral_derivative();
    //            }

    //            for(uint c = 0; c < m_basisSet.size(); c++){
    //                const BasisSet *coreC = m_basisSet.at(c);
    //                const int coreCharge = coreC->coreCharge();
    //                integrator.setCorePositionC(coreC->corePosition());

    //                if((unsigned)N == c){
    //                    differentiateWrtC = true;
    //                }else{
    //                    differentiateWrtC = false;
    //                }

    //                dhab.row(1) -= coreCharge* primitiveA.weight() * primitiveB.weight()*
    //                        integrator.nuclearAttractionIntegral_derivative(differentiateWrtA,
    //                                                                        differentiateWrtB,
    //                                                                        differentiateWrtC);
    //            }

    //        }

    //    }

    return dhab;
}


rowvec ElectronicSystem::getTwoParticleIntegralDerivative(const int a, const int b, const int c, const int d,
                                                          const int N)
{
    rowvec dQabcd = zeros<rowvec>(3);

    //    bool differentiateWrtA;
    //    bool differentiateWrtB;
    //    bool differentiateWrtC;
    //    bool differentiateWrtD;

    //    if(N == m_coreID.at(a)){
    //        differentiateWrtA = true;
    //    }else{
    //        differentiateWrtA = false;
    //    }

    //    if(N == m_coreID.at(b)){
    //        differentiateWrtB = true;
    //    }else{
    //        differentiateWrtB = false;
    //    }

    //    if(N == m_coreID.at(c)){
    //        differentiateWrtC = true;
    //    }else{
    //        differentiateWrtC = false;
    //    }

    //    if(N == m_coreID.at(d)){
    //        differentiateWrtD = true;
    //    }else{
    //        differentiateWrtD = false;
    //    }


    //    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    //    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    //    const BasisSet *coreC = m_basisSet.at(m_coreID.at(c));
    //    const BasisSet *coreD = m_basisSet.at(m_coreID.at(d));
    //    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    //    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));
    //    const ContractedGTO &contractedC = coreC->getContracted(c - m_cumSumContracted.at(m_coreID.at(c)));
    //    const ContractedGTO &contractedD = coreD->getContracted(d - m_cumSumContracted.at(m_coreID.at(d)));

    //    integrator.setCorePositionA(coreA->corePosition());
    //    integrator.setCorePositionB(coreB->corePosition());
    //    integrator.setCorePositionC(coreC->corePosition());
    //    integrator.setCorePositionD(coreD->corePosition());

    //    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
    //        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
    //        integrator.setPrimitiveA(primitiveA);

    //        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
    //            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
    //            integrator.setPrimitiveB(primitiveB);

    //            integrator.updateHermiteCoefficients(true, false, false);
    //            integrator.updateHermiteCoefficients_derivative(true, false,false);

    //            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
    //                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
    //                integrator.setPrimitiveC(primitiveC);

    //                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
    //                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
    //                    integrator.setPrimitiveD(primitiveD);

    //                    integrator.updateHermiteCoefficients(false, true, false);
    //                    integrator.updateHermiteCoefficients_derivative(false, true,false);

    //                    dQabcd += primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
    //                            * integrator.electronRepulsionIntegral_derivative(differentiateWrtA, differentiateWrtB,
    //                                                                              differentiateWrtC, differentiateWrtD);

    //                }
    //            }
    //        }
    //    }

    return dQabcd;
}
















