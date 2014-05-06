#ifndef ATOM_H
#define ATOM_H

#include "../parser/turbomoleparser.h"


namespace hf{

class Atom
{
public:
    Atom(string basisFile, const rowvec& corePosition,
         const rowvec &coreVelocity = zeros<rowvec>(3));


    const int &coreCharge() const;
    const int &angularMomentum() const;
    const int &atomType() const;
    const int &nElectrons() const;
    int nContractedGTOs() const;
    const double &coreMass() const;
    const vector<ContractedGTO> &contractedGTOs() const;


    const double &corePartialCharge() const;
    void setCorePartialCharge(const double &corePartialCharge);

    const rowvec &corePosition() const;
    void setCorePosition(const rowvec &corePosition);


    const rowvec &coreVelocity() const;
    void setCoreVelocity(const rowvec &coreVelocity);

    const string basisType() const;
    void setBasisType(const string &basisType);

    const rowvec& force() const;
    void addForce(const rowvec &coreAcceleration);

    bool frozen() const;
    void freeze(bool frozen);


private:
    int m_atomType;
    int m_nElectrons;
    int m_coreCharge;
    int m_angularMomentum;
    double m_coreMass;
    double m_corePartialCharge;
    bool m_frozen = false;
    string m_basisType;

    rowvec m_corePosition;
    rowvec m_coreVelocity;
    rowvec m_force;
    vector<ContractedGTO> m_contractedGTOs;

    void loadBasisFile(string fileName);
};

}
#endif // ATOM_H
