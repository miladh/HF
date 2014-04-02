#ifndef TURBOMOLEPARSER_H
#define TURBOMOLEPARSER_H

#include <iostream>
#include <armadillo>
#include <boost/regex.hpp>

#include"../defines.h"
#include "../contractedGTO/contractedGTO.h"

using namespace arma;
using namespace std;
using namespace boost;


namespace hf{

class TurbomoleParser
{
public:
    TurbomoleParser(string filename);

    int atomType() const;
    void setAtomType(int atomType);

    int atomCharge() const;
    void setAtomCharge(int atomCharge);

    double atomMass() const;
    void setAtomMass(double atomMass);

    int angularMomentum() const;
    void setAngularMomentum(int angularMomentum);

    vector<ContractedGTO> contractedGTOs() const;


private:
    void loadfile(string filename);
    int m_atomType;
    int m_atomCharge;
    double m_atomMass;
    int m_angularMomentum;

    vector<ContractedGTO> m_contractedGTOs;


    rowvec pow1 = {0, 0, 0};
    rowvec pow2 = {1, 0, 0};
    rowvec pow3 = {0, 1, 0};
    rowvec pow4 = {0, 0, 1};

    rowvec pow5  = {2, 0, 0};
    rowvec pow6  = {0, 2, 0};
    rowvec pow7  = {0, 0, 2};
    rowvec pow8  = {1, 1, 0};
    rowvec pow9  = {1, 0, 1};
    rowvec pow10 = {0, 1, 1};
};

}

#endif // TURBOMOLEPARSER_H
