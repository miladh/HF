#include "hfsolver.h"

using namespace hf;

HFsolver::HFsolver(ElectronicSystem *system, const int &rank, const int &nProcs):
    m_system(system),
    m_rank(rank),
    m_nProcs(nProcs),
    m_step(0),
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

    m_rank = 0;
    m_nProcs = 1;

    // MPI----------------------------------------------------------------------
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs);
#endif

    int totFunctionCalls = 0.5 * m_nBasisFunctions * (m_nBasisFunctions + 1);
    int nFunctionCallsPerProc = ceil(double(totFunctionCalls)/m_nProcs);

    m_basisIndexToProcsMap = ivec(m_nBasisFunctions);
    vector<mpiTask> taskVector;


    for (int p = 0; p < m_nBasisFunctions; p++){
        mpiTask task;
        task.nFunctionCalls = m_nBasisFunctions - p;
        task.isAvailable = true;
        task.p = p;
        taskVector.push_back(task);
    }

    for(int proc = 0; proc < m_nProcs; proc++){
        int nMyFunctionCalls = 0;

        for(mpiTask &task: taskVector){
            if(task.isAvailable
                    && nMyFunctionCalls + task.nFunctionCalls <= nFunctionCallsPerProc){
                nMyFunctionCalls += task.nFunctionCalls;
                task.isAvailable = false;

                if (m_rank == proc){
                    m_myBasisIndices.push_back(task.p);
                }
                m_basisIndexToProcsMap(task.p) = proc;
            }
        }
    }
    //---------------------------------------------------------------------------


}

void HFsolver::runSolver()
{
    double begin = MPI_Wtime();
    setupOneParticleMatrix();

    double laps = MPI_Wtime();
    setupTwoParticleMatrix();

    double end = MPI_Wtime();
    if(m_rank==0){
        cout << setprecision(3)
             << "Elapsed time on matrix setup: "<< (double(end - begin))
             << "s - " <<(double(end - laps))/(double(end - begin) +1e-10) * 100
             << "% spent on two-body term " << endl;
    }

    updateFockMatrix();

    laps = MPI_Wtime();
    advance();
    end = MPI_Wtime();

    calculateEnergy();

    if(m_rank==0){
        cout << "Elapsed time on SCF: "<< (double(end - laps)) << "s " << endl;
        cout << setprecision(14)
             << "Configuration "      << m_step
             << " - SCF iterations: " << m_iteration
             << " - Energy: "         << m_energy << endl;
        cout << "-------------------------------------------------------------------------------------"  << endl;
    }
    m_step+=1;
}



void HFsolver::setupTwoParticleMatrix()
{
    double begin = MPI_Wtime();
    for(int p: m_myBasisIndices){
        for(int r = 0; r < m_nBasisFunctions; r++){
            for(int q = p; q < m_nBasisFunctions; q++){
                for(int s = r; s < m_nBasisFunctions; s++){
                    m_Q(p,r)(q,s) = m_system->twoParticleIntegral(p,q,r,s);
                }
            }
        }
    }
    double end = MPI_Wtime();
    cout << setprecision(3)
         << "Elapsed time on two-electron integral (rank = " << m_rank << "): "
         << (double(end - begin)) << "s" << endl;


    begin = MPI_Wtime();
    for (int p = 0; p < m_nBasisFunctions; p++) {
        for (int r = 0; r < m_nBasisFunctions; r++) {
#ifdef USE_MPI
            MPI_Bcast(m_Q(p,r).memptr(), m_nBasisFunctions*m_nBasisFunctions , MPI_DOUBLE, m_basisIndexToProcsMap(p), MPI_COMM_WORLD ) ;
#endif
            for(int q = p; q < m_nBasisFunctions; q++){
                for(int s = r; s < m_nBasisFunctions; s++){
                    m_Q(q,r)(p,s) = m_Q(p,r)(q,s);
                    m_Q(p,s)(q,r) = m_Q(p,r)(q,s);
                    m_Q(q,s)(p,r) = m_Q(p,r)(q,s);
                    m_Q(r,p)(s,q) = m_Q(p,r)(q,s);
                    m_Q(s,p)(r,q) = m_Q(p,r)(q,s);
                    m_Q(r,q)(s,p) = m_Q(p,r)(q,s);
                    m_Q(s,q)(r,p) = m_Q(p,r)(q,s);
                }
            }
        }
    }

    end = MPI_Wtime();
    if(m_rank==0){
        cout <<"Communication time: "<< (double(end - begin)) <<"s" << endl;
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

const double& HFsolver::getEnergy() const
{
    return m_energy;
}

const mat& HFsolver::overlapMatrix()const
{

    return m_S;
}



