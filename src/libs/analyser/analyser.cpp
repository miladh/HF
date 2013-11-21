#include "analyser.h"

Analyser::Analyser(System *system, int nGridPoints):
    m_system(system),
    bondLength(linspace(0.5,5.0, nGridPoints)),
    m_energy(zeros(nGridPoints)),
    m_Force(zeros(nGridPoints)),
    solver(new HFsolver(system))
{
}

void Analyser::calculatePES()
{

    for(uint x = 0; x < bondLength.n_elem; x++){
        rowvec X = {bondLength(x) , 0 ,0 };
        m_system->m_basisSet.at(0)->setCorePosition(X * -0.5);
        m_system->m_basisSet.at(1)->setCorePosition(X *  0.5);

        solver->runSolver();
        m_energy(x) = solver->getEnergy();

//        cout << "[" << bondLength(x) << "," << m_energy[x] << "]," <<endl;
    }

    writeToFile(m_energy, 0);

}

void Analyser::calculateForces()
{
    for(int i = 0; i < m_energy.n_elem-1; i++){
        m_Force(i) = (m_energy(i+1) - m_energy(i) ) / (bondLength(i+1)- bondLength(i));

//        cout << "[" << bondLength(i)  << "," << m_Force(i) << "]," << endl;
    }

    writeToFile(-m_Force, 1);
}



void Analyser::writeToFile(const vec R, int n){
    stringstream outName;
    ofstream myfile;

    outName << "/home/milad/kurs/Data"<< n <<".dat";
    myfile.open(outName.str().c_str(),ios::binary);


    for(uint i=0;  i < R.n_elem; i++){
        myfile  << bondLength(i) << "  " << R(i)<< endl;
    }

    outName.str( std::string() );
    outName.clear();
    myfile.close();

}

