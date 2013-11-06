#include "system.h"
#include <basisSet/h_quadzeta.h>


System::System(int nOrbitals,int maxAngularMomentum, rowvec coreCharges, int nElectrons):
    m_h(zeros(nOrbitals,nOrbitals)),
    m_S(zeros(nOrbitals,nOrbitals)),
    m_coreCharges(coreCharges),
    m_nElectrons(nElectrons)
{

    m_Q.set_size(nOrbitals, nOrbitals);

    for(int i = 0; i < nOrbitals; i++){
        for(int j = 0; j < nOrbitals; j++){
            m_Q(i,j) = zeros(nOrbitals,nOrbitals);
        }
    }

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


mat System::getOverlapMatrix() const
{
    return m_S;
}

mat System::getOneParticleMatrix() const
{
    return m_h;
}

field<mat> System::getTwoParticleMatrix() const
{
    return m_Q;
}


void System::setupOneParticleMatrix()
{
    for(uint p = 0; p < m_coreID.size(); p++){
        for(uint q = 0; q < m_coreID.size(); q++){
            rowvec tmp = getOneParticleIntegral(p,q);
            m_S(p,q) = tmp(0);
            m_h(p,q) = tmp(1);
        }
    }

}


void System::setupTwoParticleMatrix()
{
    for(uint p = 0; p < m_coreID.size(); p++){
        for(uint r = 0; r < m_coreID.size(); r++){
            for(uint q = 0; q < m_coreID.size(); q++){
                for(uint s = 0; s < m_coreID.size(); s++){

                    m_Q(p,r)(q,s) = getTwoParticleIntegral(p,q,r,s);
                }
            }
        }
    }


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

            integrator.updateHermiteCoefficients();

            Sab += integrator.overlapIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2))
                    * primitiveA.weight() * primitiveB.weight();

            hab += integrator.kineticIntegral(powA(0), powA(1), powA(2),powB(0), powB(1), powB(2));


            for(uint c = 0; c < m_basisSet.size(); c++){
                const BasisSet *coreC = m_basisSet.at(c);
                integrator.setCorePositionC(coreC->corePosition());
                hab -= m_coreCharges(c)* integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),
                                                                              powB(0), powB(1), powB(2));

            }
            hab *= primitiveA.weight() * primitiveB.weight();

        }

    }

    rowvec tmp = {Sab, hab};
    return tmp;
}

double System::getTwoParticleIntegral(const int a, const int b, const int c, const int d)
{
    double Qabcd = 0.0;
    bool twoParticleIntegral = true;

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

            for(int k = 0; k < contractedC.getNumPrimitives(); k++){
                const PrimitiveGTO &primitiveC = contractedC.getPrimitive(k);
                const rowvec &powC = primitiveC.powers();
                integrator.setExponentC(primitiveC.exponent());

                for(int l = 0; l < contractedD.getNumPrimitives(); l++){
                    const PrimitiveGTO &primitiveD = contractedD.getPrimitive(l);
                    const rowvec &powD = primitiveD.powers();
                    integrator.setExponentD(primitiveD.exponent());

                    integrator.updateHermiteCoefficients(twoParticleIntegral);

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









//void System::setupTwoParticleMatrix()
//{

//    bool twoParticleIntegral = true;

//    for(uint A = 0; A < m_R.n_rows; A++){
//        integrator.setCorePositionA(m_R.row(A));

//        for(uint B = 0; B < m_R.n_rows; B++){
//            integrator.setCorePositionB(m_R.row(B));

//            for(uint C = 0; C < m_R.n_rows; C++){
//                integrator.setCorePositionC(m_R.row(C));

//                for(uint D = 0; D < m_R.n_rows; D++){
//                    integrator.setCorePositionD(m_R.row(D));

//                    for(uint a=0; a < m_primitives.size(); a++){
//                        integrator.setExponentA(m_primitives.at(a)->exponent());

//                        for(uint b=0; b < m_primitives.size(); b++){
//                            integrator.setExponentB(m_primitives.at(b)->exponent());

//                            for(uint c=0; c < m_primitives.size(); c++){
//                                integrator.setExponentC(m_primitives.at(c)->exponent());

//                                for(uint d=0; d < m_primitives.size(); d++){
//                                    integrator.setExponentD(m_primitives.at(d)->exponent());

//                                    integrator.updateHermiteCoefficients(twoParticleIntegral);

//                                    m_Q(a+A*4,b+B*4)(c+C*4, d+D*4) =
//                                            integrator.electronRepulsionIntegral(0,0,0,0,0,0,
//                                                                                 0,0,0,0,0,0);

//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//}



//void System::setupOneParticleMatrix()
//{
//    for(uint A = 0; A < m_R.n_rows; A++){
//        integrator.setCorePositionA(m_R.row(A));

//        for(uint B = 0; B < m_R.n_rows; B++){
//            integrator.setCorePositionB(m_R.row(B));

//            for(uint a=0; a < m_primitives.size(); a++){
//                integrator.setExponentA(m_primitives.at(a)->exponent());

//                for(uint b=0; b < m_primitives.size(); b++){
//                    integrator.setExponentB(m_primitives.at(b)->exponent());

//                    integrator.updateHermiteCoefficients();

//                    m_S(a+A*4,b+B*4) = integrator.overlapIntegral(0,0,0,0,0,0);
//                    m_h(a+A*4,b+B*4) = integrator.kineticIntegral(0,0,0,0,0,0);

//                    for(uint C = 0; C < m_R.n_rows; C++){
//                        integrator.setCorePositionC(m_R.row(C));
//                        m_h(a+A*4,b+B*4) -= integrator.nuclearAttractionIntegral(0,0,0,0,0,0);
//                    }
//                }
//            }
//        }
//    }
//}




