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
    RHF(ElectronicSystem *system);

    field<const mat *> fockMatrix();
    field<const mat *> densityMatrix() const;
    field<const mat *> expansionCoefficients() const;
    field<const vec *> orbitalEnergies() const;
    void setInitialDensity(field<mat> density);

private:
    mat m_F, m_C, m_P;
    vec m_fockEnergy;
    vector<mat> m_errors, m_fockMatrices;


protected:
    void advance();
    void solveSingle();
    void calculateEnergy();
    void updateFockMatrix();
    void DIISprocedure();
};

}

#endif // RHF_H
