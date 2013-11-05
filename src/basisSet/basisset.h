#ifndef BASISSET_H
#define BASISSET_H

#include <iostream>
#include <armadillo>

#include <src/contractedGTO/contractedGTO.h>

using namespace arma;
using namespace std;

class BasisSet
{
public:
    BasisSet();

    rowvec corePosition() const;
    void setCorePosition(const rowvec &corePosition);

    vector<ContractedGTO *> contractedGTOs() const;

    int getNumContracted() const;
protected:
    vector<ContractedGTO *> m_contractedGTOs;
    rowvec m_corePosition;

};

#endif // BASISSET_H
