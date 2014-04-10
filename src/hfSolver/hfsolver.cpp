#include "hfsolver.h"

using namespace hf;

HFsolver::HFsolver(ElectronicSystem *system):
    m_system(system),
    m_rank(0),
    m_nProcs(1),
    m_nElectrons(system->nElectrons()),
    m_nSpinUpElectrons(system->nSpinUpElectrons()),
    m_nSpinDownElectrons(system->nSpinDownElectrons()),
    m_nBasisFunctions(system->nBasisFunctions()),
    m_S(zeros(m_nBasisFunctions,m_nBasisFunctions)),
    m_h(zeros(m_nBasisFunctions,m_nBasisFunctions))

{
    m_Q.set_size(m_nBasisFunctions, m_nBasisFunctions);
    for(int i = 0; i < m_nBasisFunctions; i++){
        for(int j = 0; j < m_nBasisFunctions; j++){
            m_Q(i,j) = zeros(m_nBasisFunctions,m_nBasisFunctions);
        }
    }

    // MPI----------------------------------------------------------------------
#if USE_MPI
    m_rank   = m_world.rank();
    m_nProcs = m_world.size();
#endif

    int nPQElements = 0.5 * m_nBasisFunctions * (m_nBasisFunctions + 1);

    m_pqIndicesToProcsMap = imat(m_nBasisFunctions, m_nBasisFunctions);
    int procs = 0;
    int s = 0;
    for (int p = 0; p < m_nBasisFunctions; p++) {
        for (int q = p; q < m_nBasisFunctions; q++) {
            if (m_rank == procs){
                m_myPQIndices.push_back(pair<int, int>(p,q));
            }
            m_pqIndicesToProcsMap(p,q) = procs;
            s++;
            if(s >= BLOCK_SIZE(procs, m_nProcs, nPQElements)){
                s = 0;
                procs++;
            }
        }
    }

    //---------------------------------------------------------------------------
}

void HFsolver::runSolver()
{
    m_timer.restart();
    setupOneParticleMatrix();

    double laps = m_timer.elapsed();
    setupTwoParticleMatrix();

    double end = m_timer.elapsed();
    if(m_rank==0){
        cout << setprecision(3)
             << "Elapsed time on matrix setup: "<< end
             << "s - " <<(end - laps)/(end +1e-10) * 100
             << "% spent on two-body term " << endl;
    }

    updateFockMatrix();

    laps =  m_timer.elapsed();
    advance();
    end  = m_timer.elapsed();

    calculateEnergy();

    if(m_rank==0){
        cout << "Elapsed time on SCF: "<< end - laps << "s " << endl;
        cout << setprecision(14)
             << " - SCF iterations: " << m_iteration
             << " - Energy: "         << m_energy << endl;
        cout << "-------------------------------------------------------------------------------------"  << endl;
    }

}



void HFsolver::setupTwoParticleMatrix()
{
    double begin = m_timer.elapsed();

    int p,q;
    for(pair<int,int>pq: m_myPQIndices){
        p = pq.first;
        q = pq.second;
        for(int r = 0; r < m_nBasisFunctions; r++){
            for(int s = r; s < m_nBasisFunctions; s++){
                m_Q(p,q)(r,s) = m_system->twoParticleIntegral(p,q,r,s);
            }
        }
    }

    double end = m_timer.elapsed();
    cout << setprecision(3)
         << "Elapsed time on two-electron integral (rank = " << m_rank << "): "
         << end - begin << "s" << endl;


    begin = m_timer.elapsed();
    for (int p = 0; p < m_nBasisFunctions; p++) {
        for(int q = p; q < m_nBasisFunctions; q++){
#if USE_MPI
            boost::mpi::broadcast(m_world, m_Q(p,q).memptr(), m_Q(p,q).n_elem, m_pqIndicesToProcsMap(p,q));
#endif
            for (int r = 0; r < m_nBasisFunctions; r++) {
                for(int s = r; s < m_nBasisFunctions; s++){
                    m_Q(p,q)(s,r) = m_Q(p,q)(r,s);
                    m_Q(q,p)(r,s) = m_Q(p,q)(r,s);
                    m_Q(q,p)(s,r) = m_Q(p,q)(r,s);
                    m_Q(r,s)(p,q) = m_Q(p,q)(r,s);
                    m_Q(r,s)(q,p) = m_Q(p,q)(r,s);
                    m_Q(s,r)(p,q) = m_Q(p,q)(r,s);
                    m_Q(s,r)(q,p) = m_Q(p,q)(r,s);
                }
            }
        }
    }

    end = m_timer.elapsed();
    if(m_rank==0){
        cout <<"Communication time: "<< end - begin <<"s" << endl;
    }

}

void HFsolver::setupOneParticleMatrix()
{
    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int q = p; q < m_nBasisFunctions; q++){
            m_S(p,q) = m_system->overlapIntegral(p,q);
            m_h(p,q) = m_system->oneParticleIntegral(p,q);
        }
    }
    m_S = symmatu(m_S);
    m_h = symmatu(m_h);
}

const mat& HFsolver::normalize(mat &C, const int& HOcoeff)
{
    double norm;
    for (int i = 0; i < HOcoeff; i++){
        norm = dot(C.col(i), m_S * C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);
    }
    return C;
}

double HFsolver::computeStdDeviation(const vec& fockEnergies, const vec& fockEnergiesOld)
{
    return sum(abs(fockEnergies - fockEnergiesOld)) / fockEnergies.n_elem;
}

const double& HFsolver::energy() const
{
    return m_energy;
}

const mat& HFsolver::overlapMatrix()const
{

    return m_S;
}



