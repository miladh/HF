#include "geometricalderivative.h"

using namespace hf;

GeometricalDerivative::GeometricalDerivative(ElectronicSystem *system, HFsolver *solver):
    m_system(system),
    m_solver(solver),
    m_nBasisFunctions(system->nBasisFunctions()),
    m_rank(0),
    m_nProcs(1)
{


    // MPI----------------------------------------------------------------------
#if USE_MPI
    m_rank   = m_world.rank();
    m_nProcs = m_world.size();
#endif

    int nPQElements =  m_nBasisFunctions * m_nBasisFunctions;

    m_pqIndicesToProcsMap = imat(m_nBasisFunctions, m_nBasisFunctions);
    int procs = 0;
    int s = 0;
    for (int p = 0; p < m_nBasisFunctions; p++) {
        for (int q = 0; q < m_nBasisFunctions; q++) {
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



    if(densityMatrices.n_elem  == 1) {
        const mat& P = (*densityMatrices(0));


        int p,q;
        for(pair<int,int>pq: m_myPQIndices){
            p = pq.first;
            q = pq.second;

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

        const mat& F = (*fockMatrices(0));

        for(int c = 0; c < int(m_gradE.n_rows); c++){
            mat dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);

            int p,q;
            for(pair<int,int>pq: m_myPQIndices){
                p = pq.first;
                q = pq.second;
                mat overlapGD = m_system->overlapIntegralGD(p,q);
                dSx(p,q) += overlapGD(c,0);
                dSy(p,q) += overlapGD(c,1);
                dSz(p,q) += overlapGD(c,2);
            }

            dSx = symmatu(dSx);
            dSy = symmatu(dSy);
            dSz = symmatu(dSz);

            m_gradE(c,0) -= 0.5 * trace(P*dSx*P*F);
            m_gradE(c,1) -= 0.5 * trace(P*dSy*P*F);
            m_gradE(c,2) -= 0.5 * trace(P*dSz*P*F);

        }




    }else{
        const mat& Pu = (*densityMatrices(0));
        const mat& Pd = (*densityMatrices(1));
        const mat& Fu = (*fockMatrices(0));
        const mat& Fd = (*fockMatrices(1));


        int p,q;
        for(pair<int,int>pq: m_myPQIndices){
            p = pq.first;
            q = pq.second;

            double Pupq = Pu(p, q);
            double Pdpq = Pd(p, q);

            mat h = m_system->oneParticleIntegralGD(p,q);
            for(int c = 0; c < int(m_gradE.n_rows); c++){
                m_gradE.row(c) +=  (Pupq + Pdpq ) * h.row(c);
            }

            for (int r = 0; r < m_nBasisFunctions; r++){
                for (int s = 0; s < m_nBasisFunctions; s++){

                    double Purs = Pu(r, s);
                    double Pdrs = Pd(r, s);

                    mat J = m_system->twoParticleIntegralGD(p,q,r,s);
                    mat K = m_system->twoParticleIntegralGD(p,s,r,q);

                    for(int c = 0; c < int(m_gradE.n_rows); c++){
                        rowvec JMinusK = J.row(c) - K.row(c);
                        m_gradE.row(c) += 0.5 * JMinusK * (Pupq * Purs + Pdpq * Pdrs)
                                       +  0.5 * (Pupq * Pdrs + Pdpq * Purs) * J.row(c);
                    }
                }
            }
        }


        for(int c = 0; c < int(m_gradE.n_rows); c++){
            mat dSx = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSy = zeros(m_nBasisFunctions,m_nBasisFunctions);
            mat dSz = zeros(m_nBasisFunctions,m_nBasisFunctions);

            int p,q;
            for(pair<int,int>pq: m_myPQIndices){
                p = pq.first;
                q = pq.second;
                mat overlapGD = m_system->overlapIntegralGD(p,q);
                dSx(p,q) += overlapGD(c,0);
                dSy(p,q) += overlapGD(c,1);
                dSz(p,q) += overlapGD(c,2);
            }

            dSx = symmatu(dSx);
            dSy = symmatu(dSy);
            dSz = symmatu(dSz);

            m_gradE(c,0) -= trace(Pu*dSx*Pu*Fu) + trace(Pd*dSx*Pd*Fd);
            m_gradE(c,1) -= trace(Pu*dSy*Pu*Fu) + trace(Pd*dSy*Pd*Fd);
            m_gradE(c,2) -= trace(Pu*dSz*Pu*Fu) + trace(Pd*dSz*Pd*Fd);

        }
    }


#if USE_MPI
    m_totGradE = 0*m_gradE;
    boost::mpi::all_reduce(m_world, m_gradE.memptr(), m_gradE.n_elem, m_totGradE.memptr(), std::plus<double>());
#else
    m_totGradE = m_gradE;
#endif


    //Nuclear repulsion term
    m_totGradE  +=m_system->nuclearPotentialGD();


}


































