#include "bomd.h"

BOMD::BOMD()
{
        m_massN  = 1836.0;
        m_coreCharges = {1 , 1, 16};
        coreMass = m_massN * m_coreCharges;
        m_nElectrons = 10;
        m_basisCoreA  =  new H_QuadZeta;
        m_basisCoreB  =  new H_QuadZeta;
        m_basisCoreC  =  new O_321G;

        rowvec  A, B, C;
        double x = 1.797*cos((180-104.45) *M_PI/180.0);
        double y = 1.797*sin((180-104.45) *M_PI/180.0);
        //Water molecule
        A = { 1.797,  0.2, 0.4};
        B = {-x,  y, -0.1};
        C = { 0.0,  0.0, 0.1};

        pos = zeros(m_coreCharges.n_elem, 3);
        posNew = zeros(m_coreCharges.n_elem, 3);
        posOld = zeros(m_coreCharges.n_elem, 3);

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


void BOMD::systemConfiguration()
{
    m_nOrbitals = m_system->getTotalNumOfBasisFunc();

    m_nSteps = 400;
    m_dtn   = 1.0;
    m_dampingFactor = 0.0;


    //initialize:
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
}

void BOMD::runDynamics()
{

    for(int nStep = 0; nStep < m_nSteps; nStep++){
        writeToFile(pos,nStep);
        m_solver->runSolver();
        m_C = m_solver->getC();
        m_P = 2*m_C*m_C.t();
        m_energy = m_solver->getEnergy();


        for(uint core = 0; core < m_coreCharges.n_elem; core++){
            setupDerivativeMatrices(core);
            m_energyGradient = calculateEnergy_derivative(core);
            IntegrateCoreForwardInTime(core);
        }

        posOld = pos;
        pos = posNew;


        m_basisCoreA->setCorePosition(pos.row(0));
        m_basisCoreB->setCorePosition(pos.row(1));
        m_basisCoreC->setCorePosition(pos.row(2));


        //cout << pos << endl;
        cout << "[" << abs(pos(0,0)- pos(1,0)) << "," << m_energy << "],"<< endl;
        //cout << abs(pos(0,0)- pos(1,0)) << ","<< endl;

        //cout << "step " << nStep << " Energy: " << m_energy << endl;
    }// End of time loop nuclei

}


void BOMD::IntegrateCoreForwardInTime(int core)
{

    mat tmp = zeros(m_nOrbitals, m_nOrbitals);

    for(int i = 0; i < 3; i ++){

        for(int k = 0; k < m_nOrbitals; k++){
            for(int l = 0; l < m_nOrbitals; l++){
                tmp(k,l) = m_dS(k, l)(i);
            }
        }

        posNew(core,i) = ( 2 * pos(core,i) * coreMass(i) - posOld(core,i) *
                           ( coreMass(i) - m_dampingFactor * m_dtn * 0.5)
                           - m_dtn * m_dtn * (m_energyGradient(i) + m_energy * dot(m_C, tmp * m_C) ) )
                           / (coreMass(i) + m_dampingFactor * m_dtn * 0.5 );

    }

}

rowvec BOMD::calculateEnergy_derivative(int core)
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


void BOMD::setupDerivativeMatrices(const int core)
{

    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nOrbitals; p++){
        for(int q = 0; q < m_nOrbitals; q++){
            m_dS(p,q) = m_system->getOverlapDerivative(p,q,core);
            m_dh(p,q) = m_system->getKineticIntegralDerivative(p,q,core)
                        /*+ m_system->getAttractionIntegralDerivative(p,q,core)*/;
        }
    }

    //Set up the dQ array:
//    for(int p = 0; p < m_nOrbitals; p++){
//        for(int r = 0; r < m_nOrbitals; r++){
//            for(int q = 0; q < m_nOrbitals; q++){
//                for(int s = 0; s < m_nOrbitals; s++){

//                    m_dQ(p,r)(q,s) = m_system->getTwoParticleIntegralDerivative(p,q,r,s,core);
//                }
//            }
//        }
//    }

}


void BOMD::writeToFile(const mat R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/data"<< n <<".xyz";
    myfile.open(outName.str().c_str(),ios::binary);
    myfile << 3    << "\n";
    myfile << "Hydrogen atoms  " << "\n";

    vector<string> name;
    name.push_back("H");
    name.push_back("H");
    name.push_back("O");

    for(uint i=0;  i < R.n_rows; i++){
        myfile << name.at(i) << R.row(i);
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}


