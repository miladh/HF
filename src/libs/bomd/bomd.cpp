#include "bomd.h"

BOMD::BOMD(System *system):
    m_system(system),
    m_nCores(system->getNumOfCores()),
    m_nElectrons(system->getNumOfElectrons()),
    m_nOrbitals(system->getTotalNumOfBasisFunc()),
    m_C(ones(m_nOrbitals,m_nElectrons/2.0)),
    m_Cp(ones(m_nOrbitals,m_nElectrons/2.0)),
    m_Cm(zeros(m_nOrbitals,m_nElectrons/2.0))


{
    m_nSteps = 300;
    m_dtn   =  4.3;
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

    pos    = zeros(m_nCores, 3);
    posNew = zeros(m_nCores, 3);
    posOld = zeros(m_nCores, 3);

    for(int core = 0; core < m_nCores; core++){
        pos.row(core) = m_system->m_basisSet.at(core)->corePosition();
    }

    posOld = pos;
    m_solver = new HFsolver(m_system);
}


void BOMD::runDynamics()
{

    for(int nStep = 0; nStep < m_nSteps; nStep++){
        m_solver->runSolver();
        m_C = m_solver->getC();
        m_P = 2*m_C*m_C.t();
        m_energy = m_solver->getEnergy();

        for(int core = 0; core < m_nCores; core++){
            setupDerivativeMatrices(core);
            m_energyGradient = calculateEnergyGradient(core);
            IntegrateCoreForwardInTime(core);
        }

        writeToFile(pos,nStep);

        posOld = pos;
        pos = posNew;
        updateCorePositions();

        cout << "step " << nStep << " Energy: " << m_energy << endl;
    }

}


void BOMD::IntegrateCoreForwardInTime(int core)
{
    int coreMass = m_system->m_basisSet.at(core)->coreMass();

    mat tmp = zeros(m_nOrbitals, m_nOrbitals);

    for(int i = 0; i < 3; i ++){

        for(int k = 0; k < m_nOrbitals; k++){
            for(int l = 0; l < m_nOrbitals; l++){
                tmp(k,l) = m_dS(k, l)(i);
            }
        }

        posNew(core,i) = ( 2 * pos(core,i) * coreMass - posOld(core,i) *
                           ( coreMass - m_dampingFactor * m_dtn * 0.5)
                           - m_dtn * m_dtn * (m_energyGradient(i) + m_energy * dot(m_C, tmp * m_C) ) )
                           / (coreMass + m_dampingFactor * m_dtn * 0.5 );

    }

}

rowvec BOMD::calculateEnergyGradient(int core)
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


void BOMD::setupDerivativeMatrices(const int core)
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

void BOMD::updateCorePositions()
{

    for(int core = 0; core < m_nCores; core++){
        m_system->m_basisSet.at(core)->setCorePosition(pos.row(core));
    }

}

void BOMD::writeToFile(const mat R, int n){
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


