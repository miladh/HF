#include "system.h"
#include <iostream>

System::System(int nOrbitals, int  nNuclei ,int maxAngularMomentum):
    m_h(zeros(nOrbitals,nOrbitals)),
    m_S(zeros(nOrbitals,nOrbitals)),
    m_R(zeros(nNuclei,3))
{

    m_Q= new double***[nOrbitals];
    for (int i = 0; i < nOrbitals; ++i) {
        m_Q[i] = new double**[nOrbitals];

        for (int j = 0; j < nOrbitals; ++j){
            m_Q[i][j] = new double*[nOrbitals];

            for (int k = 0; k < nOrbitals; ++k){
                m_Q[i][j][k] = new double[nOrbitals];
            }
        }
    }
    integrator.setMaxAngularMomentum(maxAngularMomentum);

    m_R(0,0) = -0.5;
    m_R(1,0) = 0.5;

}

void System::addPrimitives(PrimitiveGTO *primitive)
{
    m_primitives.push_back(primitive);

}


mat System::getOverlapMatrix() const
{
    return m_S;
}

mat System::getOneParticleMatrix() const
{
    return m_h;
}

double**** System::getTwoParticleMatrix() const
{
    return m_Q;
}




void System::setupOneParticleMatrix()
{

    for(uint A = 0; A < m_R.n_rows; A++){
        integrator.setCorePositionA(m_R.row(A));

        for(uint B = 0; B < m_R.n_rows; B++){
            integrator.setCorePositionB(m_R.row(B));

            for(uint a=0; a < m_primitives.size(); a++){
                integrator.setExponentA(m_primitives.at(a)->exponent());

                for(uint b=0; b < m_primitives.size(); b++){
                    integrator.setExponentB(m_primitives.at(b)->exponent());
                    integrator.setupE();

                    m_S(a+A*4,b+B*4) = integrator.overlapIntegral(0,0,0,0,0,0);
                    m_h(a+A*4,b+B*4) = integrator.kineticIntegral(0,0,0,0,0,0);

                    for(uint C = 0; C < m_R.n_rows; C++){
                        integrator.setCorePositionC(m_R.row(C));
                        m_h(a+A*4,b+B*4) -= integrator.nuclearAttractionIntegral(0,0,0,0,0,0);
                    }
                }
            }
        }
    }
}

void System::setupTwoParticleMatrix()
{

    for(uint A = 0; A < m_R.n_rows; A++){
        integrator.setCorePositionA(m_R.row(A));

        for(uint C = 0; C < m_R.n_rows; C++){
            integrator.setCorePositionC(m_R.row(C));

            for(uint B = 0; B < m_R.n_rows; B++){
                integrator.setCorePositionB(m_R.row(B));

                for(uint D = 0; D < m_R.n_rows; D++){
                    integrator.setCorePositionD(m_R.row(D));

                    for(uint a=0; a < m_primitives.size(); a++){
                        integrator.setExponentA(m_primitives.at(a)->exponent());

                        for(uint c=0; c < m_primitives.size(); c++){
                            integrator.setExponentC(m_primitives.at(c)->exponent());

                            for(uint b=0; b < m_primitives.size(); b++){
                                integrator.setExponentB(m_primitives.at(b)->exponent());

                                for(uint d=0; d < m_primitives.size(); d++){
                                    integrator.setExponentD(m_primitives.at(d)->exponent());

                                    m_Q[a+A*4][c+C*4][b+B*4][d+D*4] =
                                            integrator.electronRepulsionIntegral(0,0,0,0,0,0,
                                                                                 0,0,0,0,0,0);
//                                    cout << m_Q[a+A*4][c+C*4][b+B*4][d+D*4] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}





