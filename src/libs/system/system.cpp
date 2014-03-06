#include "system.h"

using namespace hf;

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
        integrator.setPrimitiveA(primitiveA);

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            integrator.setPrimitiveB(primitiveB);

            integrator.updateHermiteCoefficients(true, false);

            Sab += integrator.overlapIntegral()
                    * primitiveA.weight() * primitiveB.weight();

            hab += primitiveA.weight() * primitiveB.weight() *
                    integrator.kineticIntegral();


            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                const int coreCharge = coreC->coreCharge();
                integrator.setCorePositionC(coreC->corePosition());
                hab -= coreCharge * primitiveA.weight() * primitiveB.weight()*
                        integrator.nuclearAttractionIntegral();

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
        integrator.setPrimitiveA(primitiveA);

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            integrator.setPrimitiveB(primitiveB);

            integrator.updateHermiteCoefficients(true, false);
            integrator.updateHermiteCoefficients_derivative(true,false);


            if(differentiateWrtA){
                dhab.row(0) += integrator.overlapIntegral_derivative()
                        * primitiveA.weight() * primitiveB.weight();

                dhab.row(1) += primitiveA.weight() * primitiveB.weight() *
                        integrator.kineticIntegral_derivative();
            }
            else if(differentiateWrtB){
                dhab.row(0) -= integrator.overlapIntegral_derivative()
                        * primitiveA.weight() * primitiveB.weight();

                dhab.row(1) -= primitiveA.weight() * primitiveB.weight() *
                        integrator.kineticIntegral_derivative();
            }

            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                const int coreCharge = coreC->coreCharge();
                integrator.setCorePositionC(coreC->corePosition());

                if((unsigned)N == c){
                    differentiateWrtC = true;
                }else{
                    differentiateWrtC = false;
                }

                dhab.row(1) -= coreCharge* primitiveA.weight() * primitiveB.weight()*
                        integrator.nuclearAttractionIntegral_derivative(differentiateWrtA,
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
        integrator.setPrimitiveA(primitiveA);

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            integrator.setPrimitiveB(primitiveB);

            integrator.updateHermiteCoefficients(true, false, false);
            integrator.updateHermiteCoefficients_derivative(true, false,false);

            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
                integrator.setPrimitiveC(primitiveC);

                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
                    integrator.setPrimitiveD(primitiveD);

                    integrator.updateHermiteCoefficients(false, true, false);
                    integrator.updateHermiteCoefficients_derivative(false, true,false);

                    dQabcd += primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
                            * integrator.electronRepulsionIntegral_derivative(differentiateWrtA, differentiateWrtB,
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
        integrator.setPrimitiveA(primitiveA);

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);
            integrator.setPrimitiveB(primitiveB);

            integrator.updateHermiteCoefficients(true, false);

            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
                integrator.setPrimitiveC(primitiveC);

                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
                    integrator.setPrimitiveD(primitiveD);

                    integrator.updateHermiteCoefficients(false, true);

                    Qabcd += primitiveA.weight() * primitiveB.weight() * primitiveC.weight() * primitiveD.weight()
                            * integrator.electronRepulsionIntegral();

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

            if((unsigned)activeCore == a){
                dVnm -= R * m_basisSet.at(a)->coreCharge() * m_basisSet.at(b)->coreCharge()/pow(dot(R,R),1.5);
            }else if((unsigned)activeCore == b){
                dVnm += R * m_basisSet.at(a)->coreCharge() * m_basisSet.at(b)->coreCharge()/pow(dot(R,R),1.5);
            }

        }
    }


    return dVnm;

}



double System::gaussianProduct(const int a, const int b, const double &x, const double &y, const double &z)
{
    double  Gab = 0.0;

    const BasisSet *coreA = m_basisSet.at(m_coreID.at(a));
    const BasisSet *coreB = m_basisSet.at(m_coreID.at(b));
    const ContractedGTO &contractedA = coreA->getContracted(a - m_cumSumContracted.at(m_coreID.at(a)));
    const ContractedGTO &contractedB = coreB->getContracted(b - m_cumSumContracted.at(m_coreID.at(b)));

    const rowvec &corePositionA = coreA->corePosition();
    const rowvec &corePositionB = coreB->corePosition();

    double Xa = x - corePositionA(0); double Xb = x - corePositionB(0);
    double Ya = y - corePositionA(1); double Yb = y - corePositionB(1);
    double Za = z - corePositionA(2); double Zb = z - corePositionB(2);

    double Ra = Xa * Xa + Ya * Ya + Za * Za;
    double Rb = Xb * Xb + Yb * Yb + Zb * Zb;

    for(int i = 0; i < contractedA.getNumPrimitives(); i++){
        const PrimitiveGTO &primitiveA = contractedA.getPrimitive(i);

        for(int j = 0; j < contractedB.getNumPrimitives(); j++){
            const PrimitiveGTO &primitiveB = contractedB.getPrimitive(j);

            Gab +=  primitiveA.weight() * primitiveB.weight()
                    * pow(Xa, primitiveA.xPower()) * pow(Xb, primitiveB.xPower())
                    * pow(Ya, primitiveA.yPower()) * pow(Yb, primitiveB.yPower())
                    * pow(Za, primitiveA.zPower()) * pow(Zb, primitiveB.zPower())
                    * exp(-primitiveA.exponent()*Ra - primitiveB.exponent()*Rb);

        }
    }

    return Gab;

}
























