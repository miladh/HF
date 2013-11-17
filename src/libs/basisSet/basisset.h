#ifndef BASISSET_H
#define BASISSET_H

#include <iostream>
#include <armadillo>
#include<includes/defines.h>
#include <contractedGTO/contractedGTO.h>

using namespace arma;
using namespace std;

class BasisSet
{
public:
    BasisSet();

    rowvec corePosition() const;
    int coreCharge() const;
    int coreMass() const;
    void setCorePosition(const rowvec &corePosition);
    void setCoreCharge(const int &coreCharge);
    void setCoreMass(const int &coreMass);

    const ContractedGTO &getContracted(const int c) const;
    int getNumContracted() const;

    virtual int getAngularMomentum() const = 0;


protected:
    vector<ContractedGTO> m_contractedGTOs;
    rowvec m_corePosition;
    int m_coreCharge, m_coreMass;
    int m_angularMomentum;

};

#endif // BASISSET_H
