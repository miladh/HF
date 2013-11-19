#include "cpmd.h"

CPMD::CPMD(System *system):
    m_system(system),
    m_nCores(system->getNumOfCores()),
    m_nElectrons(system->getNumOfElectrons()),
    m_nOrbitals(system->getTotalNumOfBasisFunc()),
    m_C(ones(m_nOrbitals,m_nElectrons/2.0)),
    m_Cp(ones(m_nOrbitals,m_nElectrons/2.0)),
    m_Cm(zeros(m_nOrbitals,m_nElectrons/2.0))
{

    m_nSteps = 300;
    m_eSteps = 50;

    m_dte    = 0.1;
    m_dtn   = 4.3;

    m_gammaE = 1.0;
    m_gammaN = 0.0;

    m_massE = 5.0;

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
            m_dh(i,j)  = zeros<rowvec>(3);
            m_dS(i,j)  = zeros<rowvec>(3);
            m_dQ(i,j).set_size(m_nOrbitals,m_nOrbitals);
        }
    }

    for(int i = 0; i < m_nOrbitals; i++ ){
        for(int j = 0; j < m_nOrbitals; j++ ){
            for(int k = 0; k < m_nOrbitals; k++ ){
                for(int l = 0; l < m_nOrbitals; l++ ){
                    m_dQ(i,j)(k,l) = zeros<rowvec>(3);
                }
            }
        }
    }

    pos    = zeros(m_nCores, 3);
    posNew = zeros(m_nCores, 3);
    posOld = zeros(m_nCores, 3);
    m_lambda = zeros<rowvec>(m_nElectrons/2.0);

    for(int core = 0; core < m_nCores; core++){
        pos.row(core) = m_system->m_basisSet.at(core)->corePosition();
    }

    posOld = pos;
    m_solver = new HFsolver(m_system);


}


void CPMD::runDynamics()
{

    for(int nStep = 0; nStep < m_nSteps; nStep++){
        m_solver->setupOneParticleMatrix();
        m_solver->setupTwoParticleMatrix();
        m_Q = m_solver->getQmatrix();
        m_h = m_solver->gethmatrix();
        m_S = m_solver->getSmatrix();

        m_C = normalize(m_C,m_S);
        m_Cm = normalize(m_Cm,m_S);

        for(uint orbital = 0; orbital < m_C.n_cols; orbital++ ){

            //Loop over time
            for(int eStep=0; eStep < m_eSteps; eStep++){
                setupFockMatrix();
                IntegrateWavefunctionForwardInTime(orbital);
                m_energy = calculateEnergy();

            }// End of time loop electrons
        }

        cout << m_lambda << endl;

        //For more then 2 electrons, lambda is different for different Cs!!!
        m_lambda *= (m_massE + m_gammaE * m_dte * 0.5 ) / ( 2 * m_dte * m_dte );

        for(int core = 0; core < m_nCores; core++){
            setupDerivativeMatrices(core);
            m_energyGradient = calculateEnergy_derivative(core);
            IntegrateCoreForwardInTime(core);
        }

        writeToFile(pos,nStep);

        posOld = pos;
        pos = posNew;
        updateCorePositions();


        cout << "step " << nStep << " Energy: " << m_energy << endl;
    }

}
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/

void CPMD::IntegrateWavefunctionForwardInTime(int orb)
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
        m_lambda(orb) = lambdaVec(0);
    }else{
        m_lambda(orb) = lambdaVec(1) ;
    }


    //Calculate C(t+h):
    m_Cp.col(orb) -= m_lambda(orb) * m_S * m_C.col(orb);


    //update C
    m_Cm.col(orb) = m_C.col(orb);
    m_C.col(orb) = m_Cp.col(orb);
}



void CPMD::IntegrateCoreForwardInTime(int core)
{
    int coreMass = m_system->m_basisSet.at(core)->coreMass();
    mat tmp = zeros(m_nOrbitals, m_nOrbitals);

    for(int i = 0; i < 3; i ++){

        for(int k = 0; k < m_nOrbitals; k++){
            for(int l = 0; l < m_nOrbitals; l++){
                tmp(k,l) = m_dS(k, l)(i);
            }
        }

        posNew(core,i) = ( 2 * pos(core,i) * coreMass - posOld(core,i) * ( coreMass - m_gammaN * m_dtn * 0.5)
                           - m_dtn * m_dtn * (m_energyGradient(i) + m_lambda(0) * dot(m_C, tmp * m_C) ) )
                         / (coreMass + m_gammaN * m_dtn * 0.5 );
    }

}

double CPMD::calculateEnergy()

{
    double Eg = 0.0;
    m_P = 2*m_C*m_C.t();
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

rowvec CPMD::calculateEnergy_derivative(int core)
{
    rowvec dE  = {0,0,0};

    for (int p = 0; p < m_nOrbitals; p++){
        for (int q = 0; q < m_nOrbitals; q++){
            dE += m_P(p, q)*m_dh(p, q);

            for (int r = 0; r < m_nOrbitals; r++){
                for (int s = 0; s < m_nOrbitals; s++){
                    dE += 0.5*m_P(p,q)*m_P(s,r)*(m_dQ(p,r)(q,s) - 0.5*m_dQ(p,r)(s,q));
                }
            }
        }
    }


    //Nuclear repulsion term
    dE  +=m_system->getNucleiPotential_derivative(core);

    return dE;
}


void CPMD::setupDerivativeMatrices(const int core)
{
    mat diffOneParticleIntegral;
    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            diffOneParticleIntegral = m_system->getOneParticleDerivative(p,q,core);
            m_dS(p,q) = diffOneParticleIntegral.row(0);
            m_dh(p,q) = diffOneParticleIntegral.row(1);
        }
    }


    //Set up the dQ array:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int r = 0; r < m_nOrbitals; r++){
            for(int q = 0; q < m_nOrbitals; q++){
                for(int s = 0; s < m_nOrbitals; s++){

                    m_dQ(p,r)(q,s) = m_system->getTwoParticleIntegralDerivative(p,q,r,s,core);
                }
            }
        }
    }


}

void CPMD::setupFockMatrix()
{

    m_P = 2*m_C*m_C.t();
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


void CPMD::updateCorePositions()
{

    for(int core = 0; core < m_nCores; core++){
        m_system->m_basisSet.at(core)->setCorePosition(pos.row(core));
    }

}

mat CPMD::normalize(mat C, mat S){

    double norm;
    for (int i = 0; i < m_nElectrons/2; i++){
        norm = dot(C.col(i), S*C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);
    }

    return C;
}


void CPMD::writeToFile(const mat R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/state"<< n <<".xyz";
    myfile.open(outName.str().c_str(),ios::binary);
    myfile << m_nCores   << "\n";
    myfile << "Hydrogen atoms  " << "\n";

    vector<string> name;
    name.push_back("H ");
    name.push_back("H ");

    if(m_nCores > 2){
         name.push_back("O ");
    }

    for(uint i=0;  i < R.n_rows; i++){
        myfile << name.at(i) << m_energy << R.row(i);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}


