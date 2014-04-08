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
    if(m_myBasisIndices.size() == 0){
         throw logic_error("Procs with no task!");
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
    for(int p: m_myBasisIndices){
        for(int r = 0; r < m_nBasisFunctions; r++){
            for(int q = p; q < m_nBasisFunctions; q++){
                for(int s = r; s < m_nBasisFunctions; s++){
                    m_Q(p,r)(q,s) = m_system->twoParticleIntegral(p,q,r,s);
                }
            }
        }
    }
    double end = m_timer.elapsed();
    cout << setprecision(3)
         << "Elapsed time on two-electron integral (rank = " << m_rank << "): "
         << end - begin << "s" << endl;


    begin = m_timer.elapsed();
    for (int p = 0; p < m_nBasisFunctions; p++) {
        for (int r = 0; r < m_nBasisFunctions; r++) {
#if USE_MPI
            boost::mpi::broadcast(m_world, m_Q(p,r).memptr(), m_Q(p,r).n_elem, m_basisIndexToProcsMap(p));
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



