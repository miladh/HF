#ifndef RHF_H
#define RHF_H

#include <armadillo>
#include <iostream>
#include "../hfSolver/hfsolver.h"

using namespace arma;
using namespace std;

namespace hf {

class RHF : public HFsolver
{
public:
    RHF(ElectronicSystem *system, const int &rank, const int &nProcs);


    const mat& getExpansionCoeff() const;
    field<mat> getFockMatrix();
    field<mat> getDensityMatrix() const;

private:
    mat m_F, m_C, m_P;
    vec m_fockEnergy;
    vector<mat> m_errors, m_fockMatrices;

    void DIISprocedure();
protected:
    void advance();
    void solveSingle();
    void calculateEnergy();
    void updateFockMatrix();
    void calculateDensity();

};

}

#endif // RHF_H
