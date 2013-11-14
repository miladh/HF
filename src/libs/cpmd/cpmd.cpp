#include "cpmd.h"

cpmd::cpmd()
{
    m_coreCharges = {1 , 1};
    m_nElectrons = 2;
    m_basisCoreA  = new H_QuadZeta;
    m_basisCoreB  = new H_QuadZeta;


    pos = zeros(m_nElectrons, 3);
    posNew = zeros(m_nElectrons, 3);
    posOld = zeros(m_nElectrons, 3);

    pos(0,0) = -0.5; pos(1,0) = 0.5;
    posOld(0,0) = -0.5; posOld(1,0) = 0.5;

    m_basisCoreA->setCorePosition(pos.row(0));
    m_basisCoreB->setCorePosition(pos.row(1));

    //NB Need higher maxL !!!!!!
    m_system = new System(m_nElectrons, m_basisCoreA->getAngularMomentum()+1, m_coreCharges);
    m_system->addBasisSet(m_basisCoreA);
    m_system->addBasisSet(m_basisCoreB);

    m_solver = new HFsolver(m_system);

    systemConfiguration();

}


void cpmd::systemConfiguration()
{
    m_nOrbitals = m_system->getTotalNumOfBasisFunc();

    m_nSteps = 100;
    m_eSteps = 400;

    m_dte    = 0.1;
    m_dtn   = 4.3;

    m_gammaE = 1.0;
    m_gammaN = 0.0;

    m_massE = 4.0;
    m_massN  = 1000 * m_massE *0.5;

    //initialize:
    m_F  = zeros(m_nOrbitals,m_nOrbitals);
    m_C  = ones(m_nOrbitals)*0.125;
    m_Cm = ones(m_nOrbitals)*0.125;
    m_Cp = zeros(m_nOrbitals);

    //   m_C  = ones(m_nOrbitals,m_nElectrons/2.0);
    //   m_Cm = ones(m_nOrbitals,m_nElectrons/2.0);
    //   m_Cp = zeros(m_nOrbitals,m_nElectrons/2.0);

    m_dh.set_size(m_nOrbitals,m_nOrbitals);
    m_dS.set_size(m_nOrbitals,m_nOrbitals);
    m_dQ.set_size(m_nOrbitals, m_nOrbitals);

    for(int i = 0; i < m_nOrbitals; i++ ){
        for(int j = 0; j < m_nOrbitals; j++ ){
            m_dh(i,j) = zeros<rowvec>(3);
            m_dS(i,j) = zeros<rowvec>(3);
            m_dQ(i,j)  = zeros(m_nOrbitals,m_nOrbitals);
        }
    }
}

/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/

void cpmd::runDynamics()
{
    double Eg;

    for(int nStep = 0; nStep < m_nSteps; nStep++){
        writeToFile(pos,nStep);
        m_solver->runSolver();
        m_Q = m_solver->getQmatrix();
        m_h = m_solver->gethmatrix();
        m_S = m_solver->getSmatrix();


        m_C = normalize(m_C,m_S);
        m_Cm = normalize(m_Cm,m_S);

        //Loop over time
        for(int eStep=0; eStep < m_eSteps+1; eStep++){
            setupFockMatrix();
            IntegrateWavefunctionForwardInTime();
            Eg = calculateEnergy(m_Q,m_h,m_C);

//            cout <<"Energy: " << setprecision(8) << Eg <<" step: " << eStep << endl;
        }// Endof time loop electrons


        m_lambda *= (m_massE + m_gammaE * m_dte * 0.5 ) / ( 2 * m_dte * m_dte );

        for(uint core = 0; core < m_coreCharges.n_elem; core++){
            setupDerivativeMatrices(core);
            m_energyGradient = calculateEnergy_derivative(m_dQ, m_dh, m_C, core);
            IntegrateCoreForwardInTime(core);
        }
        sleep(6);

        posOld = pos;
        pos = posNew;


        m_basisCoreA->setCorePosition(pos.row(0));
        m_basisCoreB->setCorePosition(pos.row(1));


               cout << pos << endl;

    }// End of time loop nuclei

}
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/


void cpmd::IntegrateCoreForwardInTime(int core)
{

    mat tmp = zeros(m_nOrbitals, m_nOrbitals);

    for(int i = 0; i < 3; i ++){

        for(int k = 0; k < m_nOrbitals; k++){
            for(int l = 0; l < m_nOrbitals; l++){
                tmp(k,l) = m_dS(k, l)(i);
            }
        }

        posNew(core,i) = ( 2 * pos(core,i) * m_massN
                           - posOld(core,i) * ( m_massN - m_gammaN * m_dtn * 0.5)
                           - m_dtn * m_dtn * (m_energyGradient(i) + m_lambda * dot(m_C, tmp * m_C) ) )
                            / (m_massN + m_gammaN * m_dtn * 0.5 );
    }


}

rowvec cpmd::calculateEnergy_derivative(const field<mat> dQ, const field<rowvec> dh, const vec C, int core)
{
    rowvec dE  = {0,0,0};
    //One-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint q=0; q < C.n_elem; q++){
            dE  += C(p)*C(q)*dh(p,q);
        }
    }
    dE = 2*dE;

    //Two-body term
    //    for(uint p=0; p < C.n_elem; p++){
    //        for(uint r=0; r < C.n_elem; r++){
    //            for(uint q=0; q< C.n_elem; q++){
    //                for(uint s=0; s < C.n_elem; s++){
    //                    dEg +=m_dQ(p,r)(q,s)*C(p)*C(q)*C(r)*C(s);
    //                }
    //            }
    //        }
    //    }


    //    Nuclear repulsion term
//        dE  +=m_system->getNucleiPotential_derivative(core);

//    rowvec R = pos.row(0) - pos.row(1);
//    dE(0) -= pow(-1.0,core)*R(0)/pow(dot(R,R),3.0/2.0);
//    dE(1) -= pow(-1.0,core)*R(1)/pow(dot(R,R),3.0/2.0);
//    dE(2) -= pow(-1.0,core)*R(2)/pow(dot(R,R),3.0/2.0);

    return dE;
}

void cpmd::setupDerivativeMatrices(const int core)
{

    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            m_dS(p,q) = m_system->getOverlapDerivative(p,q,core);
            m_dh(p,q) = /*m_system->getKineticIntegralDerivative(p,q,core)
                        +*/ m_system->getAttractionIntegralDerivative(p,q,core);
        }
    }



    mat tmp = zeros(m_nOrbitals, m_nOrbitals);
    for(int k = 0; k < m_nOrbitals; k++){
        for(int l = 0; l < m_nOrbitals; l++){
            tmp(k,l) = m_dh(k, l)(0);
        }
    }

    cout << tmp << endl;

    //    //Set up the dQ array:
    //    for(uint Rp = 0; Rp < R.n_rows; Rp++){
    //        for(uint Rr = 0; Rr < R.n_rows; Rr++){
    //            for(uint Rq = 0; Rq < R.n_rows; Rq++){
    //                for(uint Rs = 0; Rs < R.n_rows; Rs++){

    //                    for(uint p=0; p <alpha.size(); p++){
    //                        for(uint r=0; r <alpha.size(); r++){
    //                            for(uint q=0; q <alpha.size(); q++){
    //                                for(uint s=0; s <alpha.size(); s++){

    //                                    dQ[p+Rp*4][r+Rr*4][q+Rq*4][s+Rs*4] = electronInteractionIntegral_derivative(p,r,q,s,Rp,Rr,Rq,Rs,alpha,R);
    //                                }
    //                            }
    //                        }
    //                    }

    //                }
    //            }
    //        }
    //    }

}



void cpmd::IntegrateWavefunctionForwardInTime()
{

    double a, b, c;
    vec2 lambdaVec;

    m_Cp = ( 2* m_massE * m_C - (m_massE - m_gammaE * m_dte * 0.5 ) * m_Cm
             - 4 * m_dte * m_dte * m_F * m_C) / ( m_massE + m_gammaE * m_dte * 0.5);


    //Calculate lambda:
    a = dot(m_C, m_S * m_S *m_S * m_C);
    b = -2 * dot(m_S * m_Cp , m_S * m_C);
    c = dot(m_Cp, m_S * m_Cp) - 1.0;


    if(b*b - 4*a*c < 0 ){
        cerr << "Complex!" <<endl;
        cerr << a << "   "<< b << "   "<< c <<endl;
        exit (EXIT_FAILURE);

    }

    lambdaVec(0)  = (-b - sqrt(b*b - 4*a*c)) /(2*a);
    lambdaVec(1)  = (-b + sqrt(b*b - 4*a*c)) /(2*a);

    if(lambdaVec(1) < 0 && lambdaVec(0) < 0 ){
        cerr << "negative roots!!" <<endl;
        cerr << lambdaVec <<endl;
        exit (EXIT_FAILURE);

    }

    if(lambdaVec(0) >= 0.0 ){
        m_lambda = lambdaVec(0);
    }else{
        m_lambda = lambdaVec(1) ;
    }


    //Calculate C(t+h):
    m_Cp -= m_lambda * m_S * m_C;


    //update C
    m_Cm = m_C;
    m_C = m_Cp;
}


void cpmd::setupFockMatrix()
{

    //Set up the F matrix
    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){

            m_F(p,q) = m_h(p,q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    m_F(p,q) += m_Q(p,r)(q,s)*m_C(r)*m_C(s);
                }
            }
        }
    }

}

double cpmd::calculateEnergy(field<mat> Q, const mat h, const vec C)

{
    double Eg = 0.0;


    //One-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint q=0; q < C.n_elem; q++){
            Eg += C(p)*C(q)*h(p,q);
        }
    }
    Eg = 2*Eg;

    //Two-body term
    for(uint p=0; p < C.n_elem; p++){
        for(uint r=0; r < C.n_elem; r++){
            for(uint q=0; q< C.n_elem; q++){
                for(uint s=0; s < C.n_elem; s++){
                    Eg +=Q(p,r)(q,s)*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }



    //Nuclear repulsion term
    Eg +=m_system->getNucleiPotential();

    return Eg;

}




vec cpmd::normalize(vec C, mat S){
    double normFactor= 0.0;

    for(uint i= 0; i < C.n_elem; i++){
        for(uint j= 0; j < C.n_elem; j++){
            normFactor += C(i)*S(i,j)*C(j);
        }
    }
    return C/sqrt(normFactor);
}



void cpmd::writeToFile(const mat R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/data"<< n <<".xyz";
    myfile.open(outName.str().c_str(),ios::binary);
    myfile << 2    << "\n";
    myfile << "Hydrogen atoms  " << "\n";

    for(uint i=0;  i < R.n_rows; i++){
        myfile << R.row(i);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}


