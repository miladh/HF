#ifndef ATOM_H
#define ATOM_H

#include "../contractedGTO/contractedGTO.h"
#include "../parser/turbomoleparser.h"

namespace hf{

class Atom
{
public:
    Atom(string basisFile, const rowvec& corePosition);


    const int &coreCharge() const;
    const int &angularMomentum() const;
    const int &atomType() const;
    const int &nElectrons() const;
    const int nContractedGTOs() const;
    const double &coreMass() const;
    const vector<ContractedGTO> &contractedGTOs() const;


    const double &corePartialCharge() const;
    void setCorePartialCharge(const double &corePartialCharge);

    const rowvec &corePosition() const;
    void setCorePosition(const rowvec &corePosition);


private:
    int m_atomType;
    int m_nElectrons;
    int m_coreCharge;
    int m_angularMomentum;
    double m_coreMass;
    double m_corePartialCharge;

    rowvec m_corePosition;
    vector<ContractedGTO> m_contractedGTOs;

    void loadBasisFile(string fileName);
};

}
#endif // ATOM_H
