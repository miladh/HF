#include "system.h"

System::System(int nElectrons, int maxAngularMomentum):
    m_nElectrons(nElectrons)
{
    integrator.setMaxAngularMomentum(maxAngularMomentum);
    m_cumSumContracted.push_back(0);
}

int System::getTotalNumOfBasisFunc()
{
    return m_coreID.size();
}


int System::getNumOfElectrons()
{
    return m_nElectrons;
}

int System::getNumOfCores()
{
    return m_basisSet.size();
}


void System::addBasisSet(BasisSet *basisSet)
{

    int core = m_basisSet.size();
    int nGTOs = basisSet->getNumContracted();

    for(int i = 0; i <nGTOs; i++){
        m_coreID.push_back(core);
    }

    m_cumSumContracted.push_back(m_cumSumContracted.back()+ nGTOs);
    m_basisSet.push_back(basisSet);

}

rowvec System::getOneParticleIntegral(const int a, const int b)
{
    double Sab = 0;
    double hab = 0;

    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));

    integrator.setCorePositionA(coreA->corePosition());
    integrator.setCorePositionB(coreB->corePosition());


    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
        const rowvec &powA = primitiveA.powers();
        integrator.setExponentA(primitiveA.exponent());

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            const rowvec &powB = primitiveB.powers();
            integrator.setExponentB(primitiveB.exponent());

            integrator.updateHermiteCoefficients(true, false);

            Sab += integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2))
                    * primitiveA.weight() * primitiveB.weight();

            hab += primitiveA.weight() * primitiveB.weight() *
                    integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2));


            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                const int coreCharge = coreC->coreCharge();
                integrator.setCorePositionC(coreC->corePosition());
                hab -= coreCharge * primitiveA.weight() * primitiveB.weight()*
                        integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),
                                                             powB(0), powB(1), powB(2));

            }

        }

    }

    rowvec tmp = {Sab, hab};
    return tmp;
}


mat System::getOneParticleDerivative(const int a, const int b, const int N)
{
    mat dhab = zeros(2,3);
    bool differentiateWrtA;
    bool differentiateWrtB;
    bool differentiateWrtC;

    if(N == m_coreID.at(a)){
        differentiateWrtA = true;
    }else{
        differentiateWrtA = false;
    }

    if(N == m_coreID.at(b)){
        differentiateWrtB = true;
    }else{
        differentiateWrtB = false;
    }


    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));

    integrator.setCorePositionA(coreA->corePosition());
    integrator.setCorePositionB(coreB->corePosition());

    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
        const rowvec &powA = primitiveA.powers();
        integrator.setExponentA(primitiveA.exponent());

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            const rowvec &powB = primitiveB.powers();
            integrator.setExponentB(primitiveB.exponent());

            integrator.updateHermiteCoefficients(true, false);
            integrator.updateHermiteCoefficients_derivative(true,false);


            if(differentiateWrtA){
            dhab.row(0) += integrator.overlapIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2))
                        * primitiveA.weight() * primitiveB.weight();

            dhab.row(1) += primitiveA.weight() * primitiveB.weight() *
                        integrator.kineticIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2));
            }
            else if(differentiateWrtB){
                dhab.row(0) -= integrator.overlapIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2))
                            * primitiveA.weight() * primitiveB.weight();

                dhab.row(1) -= primitiveA.weight() * primitiveB.weight() *
                            integrator.kineticIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2));
            }

            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                const int coreCharge = coreC->coreCharge();
                integrator.setCorePositionC(coreC->corePosition());

                if(N == c){
                    differentiateWrtC = true;
                }else{
                    differentiateWrtC = false;
                }

                dhab.row(1) -= coreCharge* primitiveA.weight() * primitiveB.weight()*
                            integrator.nuclearAttractionIntegral_derivative(powA(0), powA(1), powA(2),
                                                                            powB(0), powB(1), powB(2),
                                                                            differentiateWrtA,
                                                                            differentiateWrtB,
                                                                            differentiateWrtC);
            }

        }

    }

    return dhab;
}


rowvec System::getTwoParticleIntegralDerivative(const int a, const int b, const int c, const int d,
                                                const int N)
{
    rowvec dQabcd = zeros<rowvec>(3);

    bool differentiateWrtA;
    bool differentiateWrtB;
    bool differentiateWrtC;
    bool differentiateWrtD;

    if(N == m_coreID.at(a)){
        differentiateWrtA = true;
    }else{
        differentiateWrtA = false;
    }

    if(N == m_coreID.at(b)){
        differentiateWrtB = true;
    }else{
        differentiateWrtB = false;
    }

    if(N == m_coreID.at(c)){
        differentiateWrtC = true;
    }else{
        differentiateWrtC = false;
    }

    if(N == m_coreID.at(d)){
        differentiateWrtD = true;
    }else{
        differentiateWrtD = false;
    }


    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    const BasisSet *coreC = m_basisSet.at(m_coreID.at(c));
    const BasisSet *coreD = m_basisSet.at(m_coreID.at(d));
    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));
    const ContractedGTO &contractedC = coreC->getContracted(c - m_cumSumContracted.at(m_coreID.at(c)));
    const ContractedGTO &contractedD = coreD->getContracted(d - m_cumSumContracted.at(m_coreID.at(d)));

    integrator.setCorePositionA(coreA->corePosition());
    integrator.setCorePositionB(coreB->corePosition());
    integrator.setCorePositionC(coreC->corePosition());
    integrator.setCorePositionD(coreD->corePosition());

    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
        const rowvec &powA = primitiveA.powers();
        integrator.setExponentA(primitiveA.exponent());

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            const rowvec &powB = primitiveB.powers();
            integrator.setExponentB(primitiveB.exponent());

            integrator.updateHermiteCoefficients(true, false);
            integrator.updateHermiteCoefficients_derivative(true, false);

            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
                const rowvec &powC = primitiveC.powers();
                integrator.setExponentC(primitiveC.exponent());

                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
                    const rowvec &powD = primitiveD.powers();
                    integrator.setExponentD(primitiveD.exponent());

                    integrator.updateHermiteCoefficients(false, true);
                    integrator.updateHermiteCoefficients_derivative(false, true);

                    dQabcd += primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
                            * integrator.electronRepulsionIntegral_derivative(powA(0), powA(1), powA(2),
                                                                              powB(0), powB(1), powB(2),
                                                                              powC(0), powC(1), powC(2),
                                                                              powD(0), powD(1), powD(2),
                                                                              differentiateWrtA, differentiateWrtB,
                                                                              differentiateWrtC, differentiateWrtD);

                }
            }
        }
    }

    return dQabcd;
}


double System::getTwoParticleIntegral(const int a, const int b, const int c, const int d)
{
    double Qabcd = 0.0;

    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    const BasisSet *coreC = m_basisSet.at(m_coreID.at(c));
    const BasisSet *coreD = m_basisSet.at(m_coreID.at(d));
    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));
    const ContractedGTO &contractedC = coreC->getContracted(c - m_cumSumContracted.at(m_coreID.at(c)));
    const ContractedGTO &contractedD = coreD->getContracted(d - m_cumSumContracted.at(m_coreID.at(d)));

    integrator.setCorePositionA(coreA->corePosition());
    integrator.setCorePositionB(coreB->corePosition());
    integrator.setCorePositionC(coreC->corePosition());
    integrator.setCorePositionD(coreD->corePosition());

    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);
        const rowvec &powA = primitiveA.powers();
        integrator.setExponentA(primitiveA.exponent());

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            const rowvec &powB = primitiveB.powers();
            integrator.setExponentB(primitiveB.exponent());

            integrator.updateHermiteCoefficients(true, false);

            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
                const rowvec &powC = primitiveC.powers();
                integrator.setExponentC(primitiveC.exponent());

                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
                    const rowvec &powD = primitiveD.powers();
                    integrator.setExponentD(primitiveD.exponent());

                    integrator.updateHermiteCoefficients(false, true);

                    Qabcd += primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
                            * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

                }
            }
        }
    }

    return Qabcd;
}


double System::getNucleiPotential()
{
    double value = 0;
    rowvec3 AB;

    for(uint a = 0; a < m_basisSet.size(); a++){
        for(uint b = a+1; b < m_basisSet.size(); b++){
            AB = m_basisSet.at(a)->corePosition() - m_basisSet.at(b)->corePosition();
            value += m_basisSet.at(a)->coreCharge() * m_basisSet.at(b)->coreCharge() /sqrt(dot(AB,AB));
        }
    }

    return value;
}

rowvec System::getNucleiPotential_derivative(int activeCore)
{
    rowvec dVnm = {0,0,0};
    rowvec3 R;



    for(uint a = 0; a < m_basisSet.size(); a++){
        for(uint b = a+1; b < m_basisSet.size(); b++){
            R = m_basisSet.at(a)->corePosition() - m_basisSet.at(b)->corePosition();

            if(activeCore == a){
                dVnm -= R * m_basisSet.at(a)->coreCharge() * m_basisSet.at(b)->coreCharge()/pow(dot(R,R),1.5);
            }else if(activeCore == b){
                dVnm += R * m_basisSet.at(a)->coreCharge() * m_basisSet.at(b)->coreCharge()/pow(dot(R,R),1.5);
            }

        }
    }


  return dVnm;

}





































//rowvec System::getTwoParticleIntegralDerivativeOfOrbitals(const int a, const int b, const int c, const int d,
//                                                const int N)
//{
//    rowvec dQabcd = zeros<rowvec>(3);

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
//        const rowvec &powA = primitiveA.powers();
//        integrator.setExponentA(primitiveA.exponent());

//        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
//            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
//            const rowvec &powB = primitiveB.powers();
//            integrator.setExponentB(primitiveB.exponent());

//            integrator.updateHermiteCoefficients(true, false);
//            integrator.updateHermiteCoefficients_derivative(true, false);

//            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
//                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
//                const rowvec &powC = primitiveC.powers();
//                integrator.setExponentC(primitiveC.exponent());

//                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
//                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
//                    const rowvec &powD = primitiveD.powers();
//                    integrator.setExponentD(primitiveD.exponent());

//                    integrator.updateHermiteCoefficients(false, true);
//                    integrator.updateHermiteCoefficients_derivative(false, true);


//                    if(differentiateWrtA){
//                        dQabcd(0) +=2 * primitiveA.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0)+1, powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powA(0)>0){
//                            dQabcd(0) -= powA(0) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0)-1, powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(1) +=2 * primitiveA.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1)+1, powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powA(1)>0){
//                            dQabcd(1) -= powA(1) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1)-1, powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(2) +=2 * primitiveA.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2)+1,powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powA(2)>0){
//                            dQabcd(2) -= powA(2) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2)-1,powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }


//                    }



//                    if(differentiateWrtB){
//                        dQabcd(0) +=2 * primitiveB.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0)+1, powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powB(0)>0){
//                            dQabcd(0) -= powB(0) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0)-1, powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(1) +=2 * primitiveB.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)+1, powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powB(1)>0){
//                            dQabcd(1) -= powB(1) * primitiveB.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)-1, powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(2) +=2 * primitiveB.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)+1,
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powB(2)>0){
//                            dQabcd(2) -= powB(2) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)-1,
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }


//                    }


//                    if(differentiateWrtC){
//                        dQabcd(0) +=2 * primitiveC.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0)+1, powC(1), powC(2),powD(0), powD(1), powD(2));

//                        if(powC(0)>0){
//                            dQabcd(0) -= powC(0) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0)-1, powC(1), powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(1) +=2 * primitiveC.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1)+1, powC(2),powD(0), powD(1), powD(2));

//                        if(powC(1)>0){
//                            dQabcd(1) -= powC(1) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1)-1, powC(2),powD(0), powD(1), powD(2));

//                        }

//                        dQabcd(2) +=2 * primitiveC.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2)+1,powD(0), powD(1), powD(2));

//                        if(powC(2)>0){
//                            dQabcd(2) -= powC(2) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2)-1,powD(0), powD(1), powD(2));

//                        }


//                    }


//                    if(differentiateWrtD){
//                        dQabcd(0) +=2 * primitiveD.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0)+1, powD(1), powD(2));

//                        if(powD(0)>0){
//                            dQabcd(0) -= powD(0) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0)-1, powD(1), powD(2));

//                        }

//                        dQabcd(1) +=2 * primitiveD.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1)+1, powD(2));

//                        if(powD(1)>0){
//                            dQabcd(1) -= powD(1) * primitiveB.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1)-1, powD(2));

//                        }

//                        dQabcd(2) +=2 * primitiveD.exponent() * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                   powC(0), powC(1), powC(2),powD(0), powD(1), powD(2)+1);

//                        if(powD(2)>0){
//                            dQabcd(2) -= powD(2) * primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
//                                    * integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2),
//                                                                       powC(0), powC(1), powC(2),powD(0), powD(1), powD(2)-1);

//                        }


//                    }




//                }
//            }
//        }
//    }

//    return dQabcd;
//}















//mat System::getOneParticleDerivativeOfOrbitals(const int a, const int b, const int N)
//{
//    mat dhab = zeros(2,3);
//    bool differentiateWrtA;
//    bool differentiateWrtB;



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
//        const rowvec &powA = primitiveA.powers();
//        integrator.setExponentA(primitiveA.exponent());

//        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
//            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
//            const rowvec &powB = primitiveB.powers();
//            integrator.setExponentB(primitiveB.exponent());

//            integrator.updateHermiteCoefficients(true, false);
//            integrator.updateHermiteCoefficients_derivative(true,false);


//            if(differentiateWrtA){
////            dhab(0,0) += 2 * primitiveA.exponent() * integrator.overlapIntegral(powA(0)+1, powA(1), powA(2),powB(0), powB(1), powB(2))
////                        * primitiveA.weight() * primitiveB.weight()
////                        - powA(0) * integrator.overlapIntegral(powA(0)-1, powA(1), powA(2),powB(0), powB(1), powB(2))
////                        * primitiveA.weight() * primitiveB.weight();

////            dhab(0,1) += 2 * primitiveA.exponent() * integrator.overlapIntegral(powA(0), powA(1)+1, powA(2),powB(0), powB(1), powB(2))
////                    * primitiveA.weight() * primitiveB.weight()
////                    - powA(1) * integrator.overlapIntegral(powA(0), powA(1)-1, powA(2),powB(0), powB(1), powB(2))
////                    * primitiveA.weight() * primitiveB.weight();

////            dhab(0,2) += 2 * primitiveA.exponent() * integrator.overlapIntegral(powA(0), powA(1), powA(2)+1,powB(0), powB(1), powB(2))
////                    * primitiveA.weight() * primitiveB.weight()
////                    - powA(2) * integrator.overlapIntegral(powA(0), powA(1), powA(2)-1,powB(0), powB(1), powB(2))
////                    * primitiveA.weight() * primitiveB.weight();



//            dhab(1,0) += 2 * primitiveA.exponent() * integrator.kineticIntegral(powA(0)+1, powA(1), powA(2),powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//            if(powA(0)> 0){
//                        dhab(1,0)-= powA(0) * integrator.kineticIntegral(powA(0)-1, powA(1), powA(2),powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//            }

//            dhab(1,1) += 2 * primitiveA.exponent() * integrator.kineticIntegral(powA(0), powA(1)+1, powA(2),powB(0), powB(1), powB(2))
//                    * primitiveA.weight() * primitiveB.weight();
//            if(powA(1)> 0){
//                    dhab(1,1)-= powA(1) * integrator.kineticIntegral(powA(0), powA(1)-1, powA(2),powB(0), powB(1), powB(2))
//                    * primitiveA.weight() * primitiveB.weight();
//            }

//            dhab(1,2) += 2 * primitiveA.exponent() * integrator.kineticIntegral(powA(0), powA(1), powA(2)+1,powB(0), powB(1), powB(2))
//                    * primitiveA.weight() * primitiveB.weight();
//            if(powA(2)> 0){
//                    dhab(1,2)-= powA(2) * integrator.kineticIntegral(powA(0), powA(1), powA(2)-1,powB(0), powB(1), powB(2))
//                    * primitiveA.weight() * primitiveB.weight();
//             }


//            for(uint c = 0; c < m_basisSet.size(); c++){
//                const BasisSet *coreC = m_basisSet.at(c);
//                const int coreCharge = coreC->coreCharge();
//                integrator.setCorePositionC(coreC->corePosition());

//                dhab(1,0) -= (2 * primitiveA.exponent() * integrator.nuclearAttractionIntegral(powA(0)+1, powA(1), powA(2),powB(0), powB(1), powB(2))
//                            * primitiveA.weight() * primitiveB.weight()) * coreCharge;

//                if(powA(0)> 0){
//                    dhab(1,0)+= powA(0) * integrator.nuclearAttractionIntegral(powA(0)-1, powA(1), powA(2),powB(0), powB(1), powB(2))
//                            * primitiveA.weight() * primitiveB.weight() * coreCharge;
//                }

//                dhab(1,1) -= 2 * primitiveA.exponent() * integrator.nuclearAttractionIntegral(powA(0), powA(1)+1, powA(2),powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight() * coreCharge;

//                if(powA(1)> 0){
//                        dhab(1,1) +=  powA(1) * integrator.nuclearAttractionIntegral(powA(0), powA(1)-1, powA(2),powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight() * coreCharge;
//                }

//                dhab(1,2) -= 2 * primitiveA.exponent() * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2)+1,powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight()* coreCharge;

//                if(powA(2)>0){
//                        dhab(1,2) +=  powA(2) * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2)-1,powB(0), powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight()* coreCharge;
//                }

//            }


//            }

//            else if(differentiateWrtB){
////                dhab(0,0) += 2 * primitiveB.exponent() * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0)+1, powB(1), powB(2))
////                        * primitiveA.weight() * primitiveB.weight()
////                        - powB(0) * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0)-1, powB(1), powB(2))
////                        * primitiveA.weight() * primitiveB.weight();

////                dhab(0,1) += 2 * primitiveB.exponent() * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)+1, powB(2))
////                        * primitiveA.weight() * primitiveB.weight()
////                        - powB(1) * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)-1, powB(2))
////                        * primitiveA.weight() * primitiveB.weight();

////                dhab(0,2) += 2 * primitiveB.exponent() * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)+1)
////                        * primitiveA.weight() * primitiveB.weight()
////                        - powB(2) * integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)-1)
////                        * primitiveA.weight() * primitiveB.weight();


//                dhab(1,0) += 2 * primitiveB.exponent() * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0)+1, powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//                if(powB(0) > 0){
//                        dhab(1,0)-= powB(0) * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0)-1, powB(1), powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//                }
//                dhab(1,1) += 2 * primitiveB.exponent() * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)+1, powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//                if(powB(1) > 0){
//                        dhab(1,1)-= powB(1) * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)-1, powB(2))
//                        * primitiveA.weight() * primitiveB.weight();
//                }
//                dhab(1,2) += 2 * primitiveB.exponent() * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)+1)
//                        * primitiveA.weight() * primitiveB.weight();
//                if(powB(2) > 0){
//                        dhab(1,2)-= powB(2) * integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)-1)
//                        * primitiveA.weight() * primitiveB.weight();
//                }


//                for(uint c = 0; c < m_basisSet.size(); c++){
//                    const BasisSet *coreC = m_basisSet.at(c);
//                    const int coreCharge = coreC->coreCharge();
//                    integrator.setCorePositionC(coreC->corePosition());

//                    dhab(1,0) -= 2 * primitiveB.exponent() * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0)+1, powB(1), powB(2))
//                            * primitiveA.weight() * primitiveB.weight()* coreCharge;

//                    if(powB(0) > 0){
//                            dhab(1,0) += powB(0) * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0)-1, powB(1), powB(2))
//                            * primitiveA.weight() * primitiveB.weight() * coreCharge;
//                    }

//                    dhab(1,1) -= 2 * primitiveB.exponent() * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)+1, powB(2))
//                            * primitiveA.weight() * primitiveB.weight()* coreCharge;
//                     if(powB(1) > 0){
//                            dhab(1,1)+= powB(1) * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1)-1, powB(2))
//                            * primitiveA.weight() * primitiveB.weight() * coreCharge;
//                    }

//                    dhab(1,2) -= 2 * primitiveB.exponent() * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)+1)
//                            * primitiveA.weight() * primitiveB.weight() * coreCharge;

//                    if(powB(2) > 0){
//                            dhab(1,2) += powB(2) * integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2)-1)
//                            * primitiveA.weight() * primitiveB.weight() * coreCharge;
//                    }


//                }



//            }

//        }

//    }

//    return dhab;
//}
