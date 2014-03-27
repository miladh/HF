#include "electronicsystem.h"

using namespace hf;

ElectronicSystem::ElectronicSystem(const int& maxAngularMomentum)
{
    integrator.setMaxAngularMomentum(maxAngularMomentum);
}

ElectronicSystem::ElectronicSystem(const int& nSpinUpElectrons, const int& nSpinDownElectrons, const int& maxAngularMomentum):
    m_nElectrons(nSpinUpElectrons + nSpinDownElectrons),
    m_nSpinUpElectrons(nSpinUpElectrons),
    m_nSpinDownElectrons(nSpinDownElectrons)

{
    integrator.setMaxAngularMomentum(maxAngularMomentum);
}

const int& ElectronicSystem::nElectrons() const
{
    return m_nElectrons;
}

const int& ElectronicSystem::nSpinUpElectrons() const
{
    return m_nSpinDownElectrons;
}

const int& ElectronicSystem::nSpinDownElectrons() const
{
    return m_nSpinDownElectrons;
}

int ElectronicSystem::nBasisFunctions()
{
    return m_basisFunctions.size();
}

int ElectronicSystem::nAtoms()
{
    return m_atoms.size();
}


void ElectronicSystem::addAtom(Atom* atom)
{

    m_atoms.push_back(atom);
    m_nElectrons += atom->nElectrons();

    for(const ContractedGTO &CGTO : atom->contractedGTOs()){
        m_basisFunctions.push_back(&CGTO);
    }
}

double ElectronicSystem::overlapIntegral(const int& p, const int& q)
{
    double Spq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    integrator.setCorePositionA(CGp->center());
    integrator.setCorePositionB(CGq->center());


    for(const PrimitiveGTO &Gp : CGp->primitivesGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitivesGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false, false);

            Spq += integrator.overlapIntegral() * Gp.weight() * Gq.weight();
        }
    }
    return Spq;
}

double ElectronicSystem::oneParticleIntegral(const int& p, const int& q)
{
    double hpq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    integrator.setCorePositionA(CGp->center());
    integrator.setCorePositionB(CGq->center());

    for(const PrimitiveGTO &Gp : CGp->primitivesGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitivesGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false,true);

            hpq += Gp.weight() * Gq.weight()* integrator.kineticIntegral();

            for(const Atom *atom : m_atoms){
                integrator.setCorePositionC(atom->corePosition());
                hpq -= atom->coreCharge() * Gp.weight() * Gq.weight() *
                        integrator.nuclearAttractionIntegral();

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
    integrator.setCorePositionA(CGp->center());
    integrator.setCorePositionB(CGq->center());
    integrator.setCorePositionC(CGr->center());
    integrator.setCorePositionD(CGs->center());

    for(const PrimitiveGTO &Gp : CGp->primitivesGTOs()) {
        integrator.setPrimitiveA(Gp);

        for(const PrimitiveGTO &Gq : CGq->primitivesGTOs()) {
            integrator.setPrimitiveB(Gq);
            integrator.updateHermiteCoefficients(true, false);

            for(const PrimitiveGTO &Gr : CGr->primitivesGTOs()) {
                integrator.setPrimitiveC(Gr);

                for(const PrimitiveGTO &Gs : CGs->primitivesGTOs()) {
                    integrator.setPrimitiveD(Gs);
                    integrator.updateHermiteCoefficients(false, true);

                    Qpqrs += Gp.weight() * Gq.weight() * Gr.weight() * Gs.weight()
                            * integrator.electronRepulsionIntegral();

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

            Vn += m_atoms.at(i)->coreCharge() * m_atoms.at(j)->coreCharge()
                    /sqrt(dot(AB,AB));
        }
    }

    return Vn;
}


double ElectronicSystem::gaussianProduct(const int& p, const int& q,
                                         const double &x, const double &y, const double &z)
{
    double  Gpq = 0.0;

    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    const rowvec &corePositionA = CGp->center();
    const rowvec &corePositionB = CGq->center();

    double Xp = x - corePositionA(0); double Xq = x - corePositionB(0);
    double Yp = y - corePositionA(1); double Yq = y - corePositionB(1);
    double Zp = z - corePositionA(2); double Zq = z - corePositionB(2);

    double Rp = Xp * Xp + Yp * Yp + Zp * Zp;
    double Rq = Xq * Xq + Yq * Yq + Zq * Zq;


    for(const PrimitiveGTO &Gp : CGp->primitivesGTOs()) {
        for(const PrimitiveGTO &Gq : CGq->primitivesGTOs()) {

            Gpq +=  Gp.weight() * Gq.weight()
                    * std::pow(Xp, Gp.xPower()) * std::pow(Xq, Gq.xPower())
                    * std::pow(Yp, Gp.yPower()) * std::pow(Yq, Gq.yPower())
                    * std::pow(Zp, Gp.zPower()) * std::pow(Zq, Gq.zPower())
                    * std::exp(-Gp.exponent()*Rp - Gq.exponent()*Rq);
        }
    }

    return Gpq;

}



void ElectronicSystem::computePartialCharge(const mat& PS)
{
    double id=0.0;


    for(Atom *atom : m_atoms){
        double partialCharge = atom->coreCharge();

        for(int j = id; j < id + atom->nContractedGTOs(); j++){
            partialCharge -= PS(j,j);
        }

        atom->setCorePartialCharge(partialCharge);
        id += atom->nContractedGTOs();
    }


    int rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if(rank == 0){
        for(Atom *atom : m_atoms){
            cout << atom->atomType() << "   " << atom->corePartialCharge() << endl;

        }

    }
}




//****************************************************************************************
//****************************************************************************************
//****************************************************************************************
//****************************************************************************************

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













