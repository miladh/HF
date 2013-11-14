#include "system.h"

System::System(int nElectrons, int maxAngularMomentum, rowvec coreCharges):
    m_coreCharges(coreCharges),
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
                integrator.setCorePositionC(coreC->corePosition());
                hab -= m_coreCharges(c)* primitiveA.weight() * primitiveB.weight()*
                        integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),
                                                             powB(0), powB(1), powB(2));

            }

        }

    }

    rowvec tmp = {Sab, hab};
    return tmp;
}


rowvec System::getOverlapDerivative(const int a, const int b, const int N)
{
    rowvec dSab = zeros<rowvec>(3);

    if(!(N == m_coreID.at(a) || N == m_coreID.at(b)) ){
        return dSab;
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

            dSab += integrator.overlapIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2))
                    * primitiveA.weight() * primitiveB.weight();

        }

    }



    if(N == m_coreID.at(a)){

        return dSab;
    }else{

        return -dSab;
    }
}


rowvec System::getKineticIntegralDerivative(const int a, const int b, const int N)
{
    rowvec dKab = zeros<rowvec>(3);

    if(!(N == m_coreID.at(a) || N == m_coreID.at(b)) ){
        return dKab;
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

            dKab += primitiveA.weight() * primitiveB.weight() *
                    integrator.kineticIntegral_derivative(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2));

        }
    }

    if(N == m_coreID.at(a)){

        return dKab;
    }else{

        return -dKab;
    }


}


rowvec System::getAttractionIntegralDerivative(const int a, const int b, const int N)
{
    rowvec dVab = zeros<rowvec>(3);

//    if(!(N == m_coreID.at(a) || N == m_coreID.at(b)) ){
//        return dVab;
//    }

    bool A;
    if(N == m_coreID.at(a)){
        A = true;
    }else{
        A = false;
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


            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                integrator.setCorePositionC(coreC->corePosition());

                dVab -= m_coreCharges(c)* primitiveA.weight() * primitiveB.weight()*
                        integrator.nuclearAttractionIntegral_derivative(powA(0), powA(1), powA(2),
                                                                        powB(0), powB(1), powB(2), A);

            }

        }
    }

    return dVab;


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
            value += m_coreCharges(a)*m_coreCharges(b)/sqrt(dot(AB,AB));
        }
    }

    return value;
}

rowvec System::getNucleiPotential_derivative(int k)
{
    double value = 0;
    rowvec dVnm = {1,1,1};
    rowvec3 R;

    for(uint j = 0; j < m_basisSet.size(); j++){
        if(k != j){
            R = m_basisSet.at(k)->corePosition() - m_basisSet.at(j)->corePosition();
            value -= 1.0/dot(R,R);
        }
    }

    for(uint i = 0; i < m_basisSet.size(); i++){
        if(k != i){
            R = m_basisSet.at(i)->corePosition() - m_basisSet.at(k)->corePosition();
            value -= 1.0/dot(R,R);
        }
    }


  return dVnm*value*0.5;

}

