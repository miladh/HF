#include "geometricalderivative.h"

using namespace hf;

GeometricalDerivative::GeometricalDerivative(ElectronicSystem *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_nBasisFunctions(system->nBasisFunctions())
{
    m_rank = 0;
    m_nProcs = 1;
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
    //---------------------------------------------------------------------------


}

const mat& GeometricalDerivative::energyGradient()
{
    calculateEnergyGradient();
    return m_totGradE;
}


void GeometricalDerivative::calculateEnergyGradient()
{
    field<const mat *> fockMatrices = m_solver->fockMatrix();
    field<const mat *> densityMatrices = m_solver->densityMatrix();
    m_gradE  = zeros(m_system->nAtoms(), 3);


    for(uint i = 0; i < densityMatrices.n_elem; i ++){
        const mat& P = (*densityMatrices(i));


        for(int p: m_myBasisIndices){
            for (int q = 0; q < m_nBasisFunctions; q++){

                mat h = m_system->oneParticleIntegralGD(p,q);
                for(int c = 0; c < int(m_gradE.n_rows); c++){
                    m_gradE.row(c) += P(p, q) * h.row(c);
                }

                for (int r = 0; r < m_nBasisFunctions; r++){
                    for (int s = 0; s < m_nBasisFunctions; s++){

                        mat J = m_system->twoParticleIntegralGD(p,q,r,s);
                        mat K = m_system->twoParticleIntegralGD(p,s,r,q);

                        for(int c = 0; c < int(m_gradE.n_rows); c++){
                            m_gradE.row(c) += 0.5*P(p,q)*P(s,r)*(J.row(c) - 0.5 * K.row(c));
                        }
                    }
                }
            }
        }
    }

#if USE_MPI
    boost::mpi::all_reduce(m_world, m_gradE.memptr(), m_gradE.n_elem, m_totGradE.memptr(), std::plus<double>());
#else
    m_totGradE = m_gradE;
#endif




    for(uint i = 0; i < densityMatrices.n_elem; i ++){
        const mat& P = (*densityMatrices(i));
        const mat& F = (*fockMatrices(i));

        for(int c = 0; c < int(m_gradE.n_rows); c++){
            mat dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);

            for(int p: m_myBasisIndices){
                for (int q = 0; q < m_nBasisFunctions; q++){
                    mat overlapGD = m_system->overlapIntegralGD(p,q);
                    dSx(p,q) += overlapGD(c,0);
                    dSy(p,q) += overlapGD(c,1);
                    dSz(p,q) += overlapGD(c,2);
                }
            }

            m_totGradE(c,0) -= 0.5 * trace(P*dSx*P*F);
            m_totGradE(c,1) -= 0.5 * trace(P*dSy*P*F);
            m_totGradE(c,2) -= 0.5 * trace(P*dSz*P*F);

        }


    }

    //    Nuclear repulsion term
    m_totGradE  +=m_system->nuclearPotentialGD();

}















//void GeometricalDerivative::calculateEnergyGradient()
//{
//    field<const mat *> fockMatrices = m_solver->fockMatrix();
//    field<const mat *> densityMatrices = m_solver->densityMatrix();
//    m_gradE  = {0,0,0};


//    for(uint i = 0; i < densityMatrices.n_elem; i ++){
//        const mat& P = (*densityMatrices(i));

//        for(int p: m_myBasisIndices){
//            for (int q = 0; q < m_nBasisFunctions; q++){
//                m_gradE += P(p, q) * m_dh(p, q);

//                for (int r = 0; r < m_nBasisFunctions; r++){
//                    for (int s = 0; s < m_nBasisFunctions; s++){
//                        m_gradE += 0.5*P(p,q)*P(s,r)*(
//                                    m_system->getTwoParticleIntegralDerivative(p,q,r,s,m_differentiationCore)
//                                    - 0.5*m_system->getTwoParticleIntegralDerivative(p,s,r,q,m_differentiationCore));

//                    }
//                }
//            }
//        }
//    }

//#if USE_MPI
//        boost::mpi::all_reduce(m_world, m_gradE.memptr(), m_gradE.n_elem, m_totGradE.memptr(), std::plus<double>());
//#endif

//    mat dSx, dSy,dSz;
//    dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
//    dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
//    dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);

//    for(int i = 0; i < m_nBasisFunctions; i++){
//        for(int j = 0; j < m_nBasisFunctions; j++){
//            dSx(i,j) = m_dS(i, j)(0);
//            dSy(i,j) = m_dS(i, j)(1);
//            dSz(i,j) = m_dS(i, j)(2);
//        }
//    }


//    for(uint i = 0; i < fockMatrices.n_elem; i++){
//        const mat& F = (*fockMatrices(i));
//        const mat& P = (*densityMatrices(i));

//        m_totGradE(0) -= 0.5 * trace(P*dSx*P*F(i));
//        m_totGradE(1) -= 0.5 * trace(P*dSy*P*F(i));
//        m_totGradE(2) -= 0.5 * trace(P*dSz*P*F(i));
//    }

//    //    Nuclear repulsion term
//    m_totGradE  +=m_system->nuclearPotentialGD(m_differentiationCore);



//}

































