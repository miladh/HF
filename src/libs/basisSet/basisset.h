#ifndef BASISSET_H
#define BASISSET_H


#include <iostream>
#include <armadillo>

#include <includes/defines.h>
#include <contractedGTO/contractedGTO.h>


using namespace arma;
using namespace std;
using namespace boost;

class BasisSet
{
public:
    BasisSet();
    BasisSet(string inFileName);

    rowvec corePosition() const;
    const int &coreCharge() const;
    const int &coreMass() const;
    void setCorePosition(const rowvec &corePosition);
    void setCoreCharge(const int &coreCharge);
    void setCoreMass(const int &coreMass);

    const ContractedGTO &getContracted(const int c) const;
    int getNumContracted() const;
    const int &getAngularMomentum() const;


private:
    vector<ContractedGTO> m_contractedGTOs;
    rowvec m_corePosition;
    int m_coreCharge, m_coreMass;
    int m_angularMomentum;

};

#endif // BASISSET_H
