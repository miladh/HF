#include "geometricalderivative.h"

using namespace hf;

GeometricalDerivative::GeometricalDerivative(System *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_nBasisFunctions(system->getTotalNumOfBasisFunc())
{
    m_dh.set_size(m_nBasisFunctions,m_nBasisFunctions);
    m_dS.set_size(m_nBasisFunctions,m_nBasisFunctions);

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


const rowvec& GeometricalDerivative::energyGradient(const int core)
{
    m_differentiationCore = core;
    setupDerivativeMatrices();
    calculateEnergyGradient();
    return m_totGradE;
}



void GeometricalDerivative::setupDerivativeMatrices()
{

    mat diffOneParticleIntegral;
    //Set up the dh and dS matrix:
    for(int p = 0; p < m_nBasisFunctions; p++){
        for(int q = 0; q < m_nBasisFunctions; q++){
            diffOneParticleIntegral = m_system->getOneParticleDerivative(p,q,m_differentiationCore);
            m_dS(p,q) = diffOneParticleIntegral.row(0);
            m_dh(p,q) = diffOneParticleIntegral.row(1);
        }
    }

}


void GeometricalDerivative::calculateEnergyGradient()
{
    const field<mat>& F = m_solver->getFockMatrix();
    const field<mat>& P = m_solver->getDensityMatrix();
    m_gradE  = {0,0,0};

    for(uint i = 0; i < F.n_elem; i ++){

        for(int p: m_myBasisIndices){
            for (int q = 0; q < m_nBasisFunctions; q++){
                m_gradE += P(i)(p, q) * m_dh(p, q);

                for (int r = 0; r < m_nBasisFunctions; r++){
                    for (int s = 0; s < m_nBasisFunctions; s++){
                        m_gradE += 0.5*P(i)(p,q)*P(i)(s,r)*(
                                    m_system->getTwoParticleIntegralDerivative(p,q,r,s,m_differentiationCore)
                                    - 0.5*m_system->getTwoParticleIntegralDerivative(p,s,r,q,m_differentiationCore));

                    }
                }
            }
        }
    }

    //MPI
    MPI_Allreduce(m_gradE.memptr(), m_totGradE.memptr(), m_gradE.n_elem, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mat dSx, dSy,dSz;
    dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
    dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
    dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);

    for(int i = 0; i < m_nBasisFunctions; i++){
        for(int j = 0; j < m_nBasisFunctions; j++){
            dSx(i,j) = m_dS(i, j)(0);
            dSy(i,j) = m_dS(i, j)(1);
            dSz(i,j) = m_dS(i, j)(2);
        }
    }


    for(uint i = 0; i < F.n_elem; i ++){
        m_totGradE(0) -= 0.5 * trace(P(i)*dSx*P(i)*F(i));
        m_totGradE(1) -= 0.5 * trace(P(i)*dSy*P(i)*F(i));
        m_totGradE(2) -= 0.5 * trace(P(i)*dSz*P(i)*F(i));
    }

    //    Nuclear repulsion term
    m_totGradE  +=m_system->getNucleiPotential_derivative(m_differentiationCore);



}

























