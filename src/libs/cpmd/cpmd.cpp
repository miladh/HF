#include "cpmd.h"

cpmd::cpmd()
{
    m_coreCharges = {1 , 1, 8};
    m_nElectrons = 10;
    m_basisCoreA  =  new H_431G;
    m_basisCoreB  =  new H_431G;
    m_basisCoreC  =  new O_431G;



    pos = zeros(m_coreCharges.n_elem, 3);
    posNew = zeros(m_coreCharges.n_elem, 3);
    posOld = zeros(m_coreCharges.n_elem, 3);

    rowvec  A, B, C;

    double x = 1.797*cos((180-104.45) *M_PI/180.0);
    double y = 1.797*sin((180-104.45) *M_PI/180.0);
    A = {1.797, 0.0, 0.0};
    B = { -x, y, 0.0};
    C = { 0.0, 0.0, 0.0};

    pos.row(0) = A;
    pos.row(1) = B;
    pos.row(2) = C;

    posOld = pos;



    m_basisCoreA->setCorePosition(pos.row(0));
    m_basisCoreB->setCorePosition(pos.row(1));
    m_basisCoreC->setCorePosition(pos.row(2));

    //NB Need higher maxL !!!!!!
    m_system = new System(m_nElectrons, m_basisCoreC->getAngularMomentum()+1, m_coreCharges);
    m_system->addBasisSet(m_basisCoreA);
    m_system->addBasisSet(m_basisCoreB);
    m_system->addBasisSet(m_basisCoreC);

    m_solver = new HFsolver(m_system);

    systemConfiguration();

}


void cpmd::systemConfiguration()
{
    m_nOrbitals = m_system->getTotalNumOfBasisFunc();

    m_nSteps = 100;
    m_eSteps = 100;

    m_dte    = 0.1;
    m_dtn   = 4.3;

    m_gammaE = 1.0;
    m_gammaN = 5.0;

    m_massE = 4.0;
    m_massN  = 1000 * m_massE *0.5;

    //initialize:
    m_F  = zeros(m_nOrbitals,m_nOrbitals);
    m_P  = zeros(m_nOrbitals,m_nOrbitals);

    m_C  = ones(m_nOrbitals,m_nElectrons/2.0);
    m_Cm = ones(m_nOrbitals,m_nElectrons/2.0);
    m_Cp = zeros(m_nOrbitals,m_nElectrons/2.0);

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

        for(uint orbital = 0; orbital < m_C.n_cols; orbital++ ){

            //Loop over time
            for(int eStep=0; eStep < m_eSteps+1; eStep++){
                setupFockMatrix();
                IntegrateWavefunctionForwardInTime(orbital);
                m_P = 2*m_C*m_C.t();
                Eg = calculateEnergy();

//                cout <<"Energy: " << setprecision(8) << Eg <<" step: " << eStep << endl;
            }// End of time loop electrons

//            cout << "orbital " << orbital << " of " << m_C.n_cols << endl;
        }


        m_lambda *= (m_massE + m_gammaE * m_dte * 0.5 ) / ( 2 * m_dte * m_dte );

        for(uint core = 0; core < m_coreCharges.n_elem; core++){
            setupDerivativeMatrices(core);
            m_energyGradient = calculateEnergy_derivative(core);
            IntegrateCoreForwardInTime(core);
        }

        posOld = pos;
        pos = posNew;


        m_basisCoreA->setCorePosition(pos.row(0));
        m_basisCoreB->setCorePosition(pos.row(1));


//               cout << pos << endl;
//        cout << "[" << abs(pos(0,0)- pos(1,0)) << "," << Eg << "],"<< endl;
//        cout << abs(pos(0,0)- pos(1,0)) << ","<< endl;

               cout << "step " << nStep << endl;
    }// End of time loop nuclei

}
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/

void cpmd::IntegrateWavefunctionForwardInTime(int orb)
{

    double a, b, c;
    vec2 lambdaVec;

    m_Cp.col(orb) = ( 2* m_massE * m_C.col(orb) - (m_massE - m_gammaE * m_dte * 0.5 ) * m_Cm.col(orb)
             - 4 * m_dte * m_dte * m_F * m_C.col(orb)) / ( m_massE + m_gammaE * m_dte * 0.5);


    //Calculate lambda:
    a = dot(m_C.col(orb), m_S * m_S *m_S * m_C.col(orb));
    b = -2 * dot(m_S * m_Cp.col(orb) , m_S * m_C.col(orb));
    c = dot(m_Cp.col(orb), m_S * m_Cp.col(orb)) - 1.0;


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
    m_Cp.col(orb) -= m_lambda * m_S * m_C.col(orb);


    //update C
    m_Cm.col(orb) = m_C.col(orb);
    m_C.col(orb) = m_Cp.col(orb);
}



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

rowvec cpmd::calculateEnergy_derivative(int core)
{
    rowvec dE  = {0,0,0};

    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){
            dE += m_P(p, q)*m_dh(p, q);

//            for (int r = 0; r < m_nOrbitals; r++){
//                for (int s = 0; s < m_nOrbitals; s++){
//                    dE += 0.5*m_P(p,q)*m_P(s,r)*(m_dQ(p,r)(q,s) - 0.5*m_dQ(p,r)(s,q));
//                }
//            }
        }
    }


    //Nuclear repulsion term
    dE  +=m_system->getNucleiPotential_derivative(core);


    return dE;
}


void cpmd::setupDerivativeMatrices(const int core)
{

    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            m_dS(p,q) = m_system->getOverlapDerivative(p,q,core);
            m_dh(p,q) = m_system->getKineticIntegralDerivative(p,q,core)
                        + m_system->getAttractionIntegralDerivative(p,q,core);
        }
    }


        //Set up the dQ array:
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

void cpmd::setupFockMatrix()
{

    //Set up the F matrix
    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){

            m_F(p,q) = m_h(p,q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    m_F(p,q) += 0.5 * m_P(s,r) * (2 * m_Q(p,r)(q,s) - m_Q(p,r)(s,q));
                }
            }
        }
    }
}

double cpmd::calculateEnergy()

{
    double Eg = 0.0;

        for (int p = 0; p < m_nOrbitals; p++){
            for (int q = 0; q < m_nOrbitals; q++){
                Eg += m_P(p, q)*m_h(p, q);

                for (int r = 0; r < m_nOrbitals; r++){
                    for (int s = 0; s < m_nOrbitals; s++){
                        Eg += 0.5*m_P(p,q)*m_P(s,r)*(m_Q(p,r)(q,s) - 0.5*m_Q(p,r)(s,q));
                    }
                }
            }
        }

    //Nuclear repulsion term
    Eg +=m_system->getNucleiPotential();

    return Eg;

}


mat cpmd::normalize(mat C, mat S){

    double norm;
    for (int i = 0; i < m_nElectrons/2; i++){
        norm = dot(C.col(i), S*C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);
    }

    return C;
}



void cpmd::writeToFile(const mat R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/dataCPMD"<< n <<".xyz";
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


