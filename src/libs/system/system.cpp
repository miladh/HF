#include "system.h"
#include<basisSet/quadzeta.h>


System::System(int nOrbitals, int  nNuclei ,int maxAngularMomentum):
    m_h(zeros(nOrbitals,nOrbitals)),
    m_S(zeros(nOrbitals,nOrbitals)),
    m_R(zeros(nNuclei,3))
{

    m_Q.set_size(nOrbitals, nOrbitals);

    for(int i = 0; i < nOrbitals; i++){
        for(int j = 0; j < nOrbitals; j++){
            m_Q(i,j) = zeros(nOrbitals,nOrbitals);
        }
    }

    integrator.setMaxAngularMomentum(maxAngularMomentum);
    m_R(0,0) = -0.5;
    m_R(1,0) = 0.5;

    m_cumSumContracted.push_back(0);

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
    for(uint a = 0; a < m_coreID.size(); a++){

        integrator.setCorePositionA(m_basisSet.at(m_coreID.at(a))->corePosition());
        vector<ContractedGTO *> cGTOsA = m_basisSet.at(m_coreID.at(a))->contractedGTOs();
        ContractedGTO *cGTOA = cGTOsA.at(a - m_cumSumContracted.at(m_coreID.at(a)) );

        for(uint b = 0; b < m_coreID.size(); b++){

            integrator.setCorePositionB(m_basisSet.at(m_coreID.at(b))->corePosition());
            vector<ContractedGTO *> cGTOsB = m_basisSet.at(m_coreID.at(b))->contractedGTOs();
            ContractedGTO *cGTOB = cGTOsB.at(b - m_cumSumContracted.at(m_coreID.at(b)) );

            for(int i = 0; i < cGTOA->getNumPrimitives(); i++){
                vector<PrimitiveGTO *>psA = cGTOA->primitives();
                PrimitiveGTO* pA = psA.at(i);
                integrator.setExponentA(pA->exponent());

                for(int j = 0; j < cGTOB->getNumPrimitives(); j++){
                    vector<PrimitiveGTO *>psB = cGTOB->primitives();
                    PrimitiveGTO* pB = psB.at(j);
                    integrator.setExponentB(pB->exponent());

                    rowvec powA = pA->powers();
                    rowvec powB = pB->powers();

                    integrator.updateHermiteCoefficients();

                    m_S(a, b) += integrator.overlapIntegral(powA(0), powA(1), powA(2),
                                                            powB(0), powB(1), powB(2));

                    m_h(a, b) += integrator.kineticIntegral(powA(0), powA(1), powA(2),
                                                            powB(0), powB(1), powB(2));


                    for(uint c = 0; c < m_basisSet.size(); c++){
                        integrator.setCorePositionC(m_basisSet.at(c)->corePosition());
                        m_h(a,b) -= integrator.nuclearAttractionIntegral(powA(0), powA(1), powA(2),
                                                                         powB(0), powB(1), powB(2));
                    }

                }

            }

        }

    }
}




void System::setupTwoParticleMatrix()
{

    bool twoParticleIntegral = true;
    for(uint a = 0; a < m_coreID.size(); a++){

        integrator.setCorePositionA(m_basisSet.at(m_coreID.at(a))->corePosition());
        vector<ContractedGTO *> cGTOsA = m_basisSet.at(m_coreID.at(a))->contractedGTOs();
        ContractedGTO *cGTOA = cGTOsA.at(a - m_cumSumContracted.at(m_coreID.at(a)) );

        for(uint b = 0; b < m_coreID.size(); b++){

            integrator.setCorePositionB(m_basisSet.at(m_coreID.at(b))->corePosition());
            vector<ContractedGTO *> cGTOsB = m_basisSet.at(m_coreID.at(b))->contractedGTOs();
            ContractedGTO *cGTOB = cGTOsB.at(b - m_cumSumContracted.at(m_coreID.at(b)) );


            for(uint c = 0; c < m_coreID.size(); c++){

                integrator.setCorePositionC(m_basisSet.at(m_coreID.at(c))->corePosition());
                vector<ContractedGTO *> cGTOsC = m_basisSet.at(m_coreID.at(c))->contractedGTOs();
                ContractedGTO *cGTOC = cGTOsC.at(c - m_cumSumContracted.at(m_coreID.at(c)) );

                for(uint d = 0; d < m_coreID.size(); d++){

                    integrator.setCorePositionD(m_basisSet.at(m_coreID.at(d))->corePosition());
                    vector<ContractedGTO *> cGTOsD = m_basisSet.at(m_coreID.at(d))->contractedGTOs();
                    ContractedGTO *cGTOD = cGTOsD.at(d - m_cumSumContracted.at(m_coreID.at(d)) );

                    for(int i = 0; i < cGTOA->getNumPrimitives(); i++){
                        vector<PrimitiveGTO *>psA = cGTOA->primitives();
                        PrimitiveGTO* pA = psA.at(i);
                        integrator.setExponentA(pA->exponent());

                        for(int j = 0; j < cGTOB->getNumPrimitives(); j++){
                            vector<PrimitiveGTO *>psB = cGTOB->primitives();
                            PrimitiveGTO* pB = psB.at(j);
                            integrator.setExponentB(pB->exponent());


                            for(int k = 0; k < cGTOC->getNumPrimitives(); k++){
                                vector<PrimitiveGTO *>psC = cGTOC->primitives();
                                PrimitiveGTO* pC = psC.at(k);
                                integrator.setExponentC(pC->exponent());

                                for(int l = 0; l < cGTOD->getNumPrimitives(); l++){
                                    vector<PrimitiveGTO *>psD = cGTOD->primitives();
                                    PrimitiveGTO* pD = psD.at(l);
                                    integrator.setExponentD(pD->exponent());

                                    rowvec powA = pA->powers();
                                    rowvec powB = pB->powers();
                                    rowvec powC = pC->powers();
                                    rowvec powD = pD->powers();

                                    integrator.updateHermiteCoefficients(twoParticleIntegral);

                                    m_Q(a,b)(c, d) =
                                            integrator.electronRepulsionIntegral(powA(0), powA(1), powA(2),
                                                                                 powB(0), powB(1), powB(2),
                                                                                 powC(0), powC(1), powC(2),
                                                                                 powD(0), powD(1), powD(2));

                                }
                            }
                        }
                    }
                }
            }
        }
    }
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




