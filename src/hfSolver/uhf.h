#ifndef UHF_H
#define UHF_H

#include <armadillo>
#include <iostream>
#include "../hfSolver/hfsolver.h"

using namespace arma;
using namespace std;

namespace hf {

class UHF : public HFsolver
{
public:
    UHF(ElectronicSystem *system);
    field<const mat *> fockMatrix();
    field<const mat *> densityMatrix() const;
    field<const mat *> expansionCoefficients() const;


private:
    mat m_Fu, m_Cu, m_Pu;
    mat m_Fd, m_Cd, m_Pd;
    vec m_fockEnergyU, m_fockEnergyD;
    vector<mat> m_errorsU, m_fockMatricesU;
    vector<mat> m_errorsD, m_fockMatricesD;

    void DIISprocedure();
protected:
    void advance();
    void solveSingle();
    void calculateEnergy();
    void updateFockMatrix();


};

}
#endif // UHF_H
