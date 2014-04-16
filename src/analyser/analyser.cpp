#include "analyser.h"
#include "../hfSolver/rhf.h"
#include "../hfSolver/uhf.h"

using namespace hf;


Analyser::Analyser(ElectronicSystem* system, HFsolver* solver, bool saveResults):
    m_system(system),
    m_solver(solver),
    m_integrator(new Integrator(m_system->maxAngularMomentum())),
    m_basisFunctions(system->basisFunctions()),
    m_nBasisFunctions(m_basisFunctions.size()),
    m_rank(0),
    m_nProcs(1)
{

#if USE_MPI
    m_rank = m_world.rank();
    m_nProcs = m_world.size();
#endif

}

Analyser::Analyser(const Config* cfg, ElectronicSystem* system, HFsolver* solver):
    m_cfg(cfg),
    m_system(system),
    m_solver(solver),
    m_integrator(new Integrator(m_system->maxAngularMomentum())),
    m_basisFunctions(system->basisFunctions()),
    m_nBasisFunctions(m_basisFunctions.size()),
    m_rank(0),
    m_nProcs(1)
{

#if USE_MPI
    m_rank = m_world.rank();
    m_nProcs = m_world.size();
#endif
    if(m_rank == 0){
        cout << "Starting analysis......" << endl;
    }
}


void Analyser::runAnalysis()
{
    const Setting & root = m_cfg->getRoot();

    if(int(root["analysisSettings"]["saveResults"])){
        string outputFilePath = root["analysisSettings"]["outputFilePath"];
        m_outputManager = new OutputManager(m_system->nAtoms(), outputFilePath);

        if(int(root["analysisSettings"]["saveEnergies"])){
            saveEnergies();
            if(m_rank == 0){
                cout << " - Energies saved!" << endl;
            }
        }
        if(int(root["analysisSettings"]["dipoleMoment"])){
            computeDipoleMoment();
            if(m_rank == 0){
                cout << " - Dipole moment computed!" << endl;
            }
        }
        if(int(root["analysisSettings"]["atomicPartialCharge"])){
            computeAtomicPartialCharge();
            if(m_rank == 0){
                cout << " - Partial charge computed!" << endl;
            }
        }

        saveResults();
    }
    else if(m_rank == 0){
        cout << " - No output file will be created!" << endl;
    }


    if(int(root["analysisSettings"]["chargeDensity"])){
        computeChargeDensity();
        if(m_rank == 0){
            cout << " - Charge density computed!" << endl;
        }
    }
    if(int(root["analysisSettings"]["electrostaticPotential"])){
        computeElectrostaticPotential();
        if(m_rank == 0){
            cout << " - Electrostatic potential computed!" << endl;
        }
    }

}

void Analyser::saveResults()
{
    m_outputManager->saveAtoms(m_system->atoms());
    m_outputManager->closeOutput();
}

void Analyser::saveEnergies()
{

    field<const vec *> energies= m_solver->orbitalEnergies();
    mat orbitalEnergies = zeros(energies.n_elem, energies.at(0)->n_elem);

    for(int i = 0; i < signed(energies.n_elem); i++){
        orbitalEnergies.row(i) = (*energies(i)).t();
    }

    m_outputManager->saveEnergy(m_solver->energy(), orbitalEnergies);

}

void Analyser::computeAtomicPartialCharge()
{
    const mat& S = m_solver->overlapMatrix();
    field<const mat *> densityMatrices = m_solver->densityMatrix();
    field<mat> PS(densityMatrices.n_elem, 1);

    for(uint i = 0; i < densityMatrices.n_elem; i ++){
        const mat& P = (*densityMatrices(i));
        PS(i) = P * S;
    }

    double id=0.0;
    for(Atom *atom : m_system->atoms()){
        double partialCharge = atom->coreCharge();

        for(int j = id; j < id + atom->nContractedGTOs(); j++){
            for(const mat& PS_i : PS)
                partialCharge -= PS_i(j,j);
        }

        atom->setCorePartialCharge(partialCharge);
        id += atom->nContractedGTOs();
    }
}


void Analyser::computeDipoleMoment()
{
    rowvec D ={0,0,0};
    field<const mat *> densityMatrices = m_solver->densityMatrix();
    const mat& P = (*densityMatrices(0));

    for(int p = 0; p < m_system->nBasisFunctions(); p++){
        for(int q = 0; q < m_system->nBasisFunctions(); q++){
            const ContractedGTO *CGp = m_basisFunctions.at(p);
            const ContractedGTO *CGq = m_basisFunctions.at(q);

            for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
                m_integrator->setPrimitiveA(Ga);

                for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
                    m_integrator->setPrimitiveB(Gb);
                    m_integrator->updateOverlapHermiteCoefficients();

                    D -= P(p,q) * m_integrator->dipoleIntegral();

                }
            }
        }
    }
    for(const Atom *atom : m_system->atoms()){
        m_integrator->setNuclearSourceCharge(atom->corePosition());
        D +=  atom->coreCharge() * atom->corePosition();

    }

    m_outputManager->saveDipoleMoment(sqrt(dot(D,D)));
}


void Analyser::computeElectrostaticPotential()
{
    vec x = linspace(-10, 10, 40);
    vec y = linspace(-10, 10, 40);
    vec z = linspace(-10, 10, 40);

    cube densityCube = zeros(y.n_elem, x.n_elem, z.n_elem);
    field<const mat *> densityMatrices = m_solver->densityMatrix();



    //MPI---------------------------------------------------------------------
    vector<int> myGridPoints;
    int node = 0;
    int s = 0;
    for (int i = 0; i < int(x.n_elem); i++) {
        if (m_rank == node){
            myGridPoints.push_back(i);
        }
        s++;
        if(s >= signed(BLOCK_SIZE(node, m_nProcs, x.n_elem))){
            s = 0;
            node++;
        }
    }

    //-----------------------------------------------------------------------
    for(uint d = 0; d < densityMatrices.n_elem; d++){
        const mat& P = (*densityMatrices(d));

        for(const int i : myGridPoints) {
            cout << m_rank << "   " << x(i) << endl;
            for(int j = 0; j < int(y.n_elem); j++) {
                for(int k = 0; k < int(z.n_elem); k++) {


                    //----------------------------------------------------------------//

                    rowvec C = {x(i), y(j), z(k)};
                    for(Atom* atom : m_system->atoms()){
                        double rAP = sqrt(dot(atom->corePosition() - C, atom->corePosition() - C));
                        densityCube(j,i,k) += atom->coreCharge() / rAP;
                    }

                    for(int p = 0; p < m_nBasisFunctions; p++){
                        densityCube(j,i,k) -= P(p,p) * electronicPotential(p,p,C);

                        for(int q = p+1; q < m_nBasisFunctions; q++){
                            densityCube(j,i,k) -= 2.0 *  P(p,q) * electronicPotential(p,q,C);

                        }
                    }

                    //------------------------------------------------------------------//


                }
            }
        }
    }

    writeDensityToFile(densityCube, x.min(),x.max(),y.min(),y.max(),z.min(),z.max());

}



double Analyser::electronicPotential(const int& p, const int& q, const rowvec& C)
{
    double Vpq = 0;
    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    m_integrator->setNuclearSourceCharge(C);

    for(const PrimitiveGTO &Ga : CGp->primitiveGTOs()) {
        m_integrator->setPrimitiveA(Ga);

        for(const PrimitiveGTO &Gb : CGq->primitiveGTOs()) {
            m_integrator->setPrimitiveB(Gb);
//            m_integrator->updateKineticHermiteCoefficients();
            m_integrator->updateOverlapHermiteCoefficients();

            Vpq += m_integrator->nuclearAttractionIntegral();

        }
    }

    return Vpq;
}


void Analyser::computeChargeDensity()
{
    vec x = linspace(-10, 10, 50);
    vec y = linspace(-10, 10, 50);
    vec z = linspace(-10, 10, 50);
    double dr = (x(1) - x(0)) * (y(1) - y(0)) * (z(1) - z(0));

    cube densityCube = zeros(y.n_elem, x.n_elem, z.n_elem);
    field<const mat *> densityMatrices = m_solver->densityMatrix();

    //MPI---------------------------------------------------------------------
    vector<int> myGridPoints;
    int node = 0;
    int s = 0;
    for (int i = 0; i < int(x.n_elem); i++) {
        if (m_rank == node){
            myGridPoints.push_back(i);
        }
        s++;
        if(s >= signed(BLOCK_SIZE(node, m_nProcs, x.n_elem))){
            s = 0;
            node++;
        }
    }
    //-----------------------------------------------------------------------

    double sumDensity =0;
    for(uint p = 0; p < densityMatrices.n_elem; p++){
        const mat& P = (*densityMatrices(p));

        for(const int i : myGridPoints) {
            for(int j = 0; j < int(y.n_elem); j++) {
                for(int k = 0; k < int(z.n_elem); k++) {

                    for(int p = 0; p < m_system->nBasisFunctions(); p++){
                        double innerProduct = gaussianProduct(p, p, x(i), y(j), z(k));
                        sumDensity += P(p,p) * innerProduct * dr;
                        densityCube(j,i,k) += P(p,p) * innerProduct;

                        for(int q = p+1; q < m_system->nBasisFunctions(); q++){
                            innerProduct = gaussianProduct(p, q, x(i), y(j), z(k));
                            sumDensity += 2.0 * P(p,q) * innerProduct * dr;
                            densityCube(j,i,k) += 2.0 * P(p,q) * innerProduct ;

                        }
                    }

                }
            }
        }
    }

    double nElectrons = sumDensity;
#if USE_MPI
    if(m_rank == 0){
        boost::mpi::reduce(m_world, sumDensity, nElectrons, std::plus<double>(), 0);
    }else{
        boost::mpi::reduce(m_world, sumDensity, std::plus<double>(), 0);
    }
#endif
    if(m_rank == 0){
        cout << "Integrated number of electrons: " << nElectrons << endl;
    }

    writeDensityToFile(densityCube, x.min(),x.max(),y.min(),y.max(),z.min(),z.max());

}


double Analyser::gaussianProduct(const int& p, const int& q,
                                 const double &x, const double &y, const double &z)
{
    double  Gpq = 0.0;

    const ContractedGTO *CGp = m_basisFunctions.at(p);
    const ContractedGTO *CGq = m_basisFunctions.at(q);
    const rowvec &corePositionA = CGp->center();
    const rowvec &corePositionB = CGq->center();

    double Xp = x - corePositionA(0); double Xq = x - corePositionB(0);
    double Yp = y - corePositionA(1); double Yq = y - corePositionB(1);
    double Zp = z - corePositionA(2); double Zq = z - corePositionB(2);

    double Rp = Xp * Xp + Yp * Yp + Zp * Zp;
    double Rq = Xq * Xq + Yq * Yq + Zq * Zq;

    for(const PrimitiveGTO &Gp : CGp->primitiveGTOs()) {
        for(const PrimitiveGTO &Gq : CGq->primitiveGTOs()) {

            Gpq +=  Gp.weight() * Gq.weight()
                    * std::pow(Xp, Gp.xPower()) * std::pow(Xq, Gq.xPower())
                    * std::pow(Yp, Gp.yPower()) * std::pow(Yq, Gq.yPower())
                    * std::pow(Zp, Gp.zPower()) * std::pow(Zq, Gq.zPower())
                    * std::exp(-Gp.exponent()*Rp - Gq.exponent()*Rq);
        }
    }

    return Gpq;

}




void Analyser::writeDensityToFile(const cube& density,
                                  const double &xMin, const double &xMax,
                                  const double &yMin, const double &yMax,
                                  const double &zMin, const double &zMax)
{

    stringstream cubeFileName;
    cubeFileName<<"/home/milad/kurs/qmd/density/id"<< m_rank <<"_cubeFile" << setw(4) << setfill('0') << 0 <<".cube";
    ofstream cubeFile(cubeFileName.str(), ios::out | ios::binary);

    //Header
    double nCores = m_system->nAtoms();
    vec origo = {0,0,0};
    vec xLimits = {xMin, xMax};
    vec yLimits = {yMin, yMax};
    vec zLimits = {zMin, zMax};
    double nX = density.n_rows;
    double nY = density.n_cols;
    double nZ = density.n_slices;

    cubeFile.write(reinterpret_cast<const char*>(&nCores), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(1)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&origo(2)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nX), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&xLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&xLimits(1)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nY), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&yLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&yLimits(1)), sizeof(double));

    cubeFile.write(reinterpret_cast<const char*>(&nZ), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&zLimits(0)), sizeof(double));
    cubeFile.write(reinterpret_cast<const char*>(&zLimits(1)), sizeof(double));

    for(const Atom* atom : m_system->atoms()){
        rowvec3 R = atom->corePosition();
        double atomType = atom->atomType();

        cubeFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&atomType), sizeof(double));

        // Write the x, y and z-components
        cubeFile.write(reinterpret_cast<const char*>(&R(0)), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&R(1)), sizeof(double));
        cubeFile.write(reinterpret_cast<const char*>(&R(2)), sizeof(double));
    }

    //density data
    for(uint k = 0; k < density.n_slices; k++) {
        for(uint i = 0; i < density.n_rows; i++) {
            for(uint j = 0; j < density.n_cols; j++) {
                double p = density(i,j,k);
                cubeFile.write(reinterpret_cast<const char*>(&p), sizeof(double));
            }
        }
    }

    cubeFile.close();
}


