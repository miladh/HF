#include "basisset.h"


BasisSet::BasisSet()
{
}

BasisSet::BasisSet(string inFileName)
{
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

    string line, stringToSearch;
    fstream file;
    file.open(inFileName);


    //Read turbomole file
    if (file.is_open()){
        while (getline(file, line)){
            stringToSearch += line +"\n";
        }
        file.close();
    } else {
        cout << "Error: Could not open file "<< inFileName << endl;
        exit(EXIT_FAILURE);
    }

    //Search for groups of contracted GTOs
    regex searchCGTO("([0-9]\\s+[spd])((\\s+-?[0-9]+\\.[0-9]+)+)");
    sregex_iterator CGTOs(stringToSearch.begin(), stringToSearch.end(), searchCGTO);
    sregex_iterator endCGTOs;

    for(; CGTOs!=endCGTOs; CGTOs++){
        string subString = CGTOs->str(2).c_str();
        //        cout << CGTOs->str(2).c_str() << endl;

        //Search for groups of (exponent and coefficient) for each primitive
        regex searchPGTO("(-?[0-9]+\\.[0-9]+)\\s+(-?[0-9]+\\.[0-9]+)\\s*");
        sregex_iterator PGTOs(subString.begin(), subString.end(), searchPGTO);
        sregex_iterator endPGTOs;

        // Check type of orbital (s, p, d)
        string orbitalType = CGTOs->str(1).c_str();
        regex searchs("[0-9]\\s+s");
        regex searchp("[0-9]\\s+p");
        regex searchd("[0-9]\\s+d");


        if(regex_match(orbitalType, searchs)){
            m_angularMomentum = 0;
            ContractedGTO contractedGTO;

            for(; PGTOs!=endPGTOs; PGTOs++){
                double exponent = atof(PGTOs->str(1).c_str());
                double weight   = atof(PGTOs->str(2).c_str());
                weight = pow(2*exponent / M_PI, 0.75) * weight;

                PrimitiveGTO primitiveGTO(exponent,weight,pow1);
                contractedGTO.addPrimitive(primitiveGTO);
            }

            m_contractedGTOs.push_back(contractedGTO);
        }
        else if(regex_match(orbitalType, searchp)){
            m_angularMomentum = 1;
            ContractedGTO contractedGTOx, contractedGTOy, contractedGTOz ;

            for(; PGTOs!=endPGTOs; PGTOs++){
                double exponent = atof(PGTOs->str(1).c_str());
                double weight   = atof(PGTOs->str(2).c_str());
                weight = pow(2*exponent / M_PI, 0.75) * 2 * sqrt(exponent) * weight;

                PrimitiveGTO primitiveGTOx(exponent,weight,pow2);
                PrimitiveGTO primitiveGTOy(exponent,weight,pow3);
                PrimitiveGTO primitiveGTOz(exponent,weight,pow4);

                contractedGTOx.addPrimitive(primitiveGTOx);
                contractedGTOy.addPrimitive(primitiveGTOy);
                contractedGTOz.addPrimitive(primitiveGTOz);
            }

            m_contractedGTOs.push_back(contractedGTOx);
            m_contractedGTOs.push_back(contractedGTOy);
            m_contractedGTOs.push_back(contractedGTOz);
        }
        else if(regex_match(orbitalType, searchd)){
            m_angularMomentum = 2;
            ContractedGTO contractedGTOxx, contractedGTOyy, contractedGTOzz,
                          contractedGTOxy, contractedGTOxz, contractedGTOyz ;

            for(; PGTOs!=endPGTOs; PGTOs++){
                double exponent = atof(PGTOs->str(1).c_str());
                double weight   = atof(PGTOs->str(2).c_str());
                weight = pow(2*exponent / M_PI, 0.75) * 2 * sqrt(exponent) * weight;

                PrimitiveGTO primitiveGTOxx(exponent,weight,pow5);
                PrimitiveGTO primitiveGTOyy(exponent,weight,pow6);
                PrimitiveGTO primitiveGTOzz(exponent,weight,pow7);
                PrimitiveGTO primitiveGTOxy(exponent,weight,pow8);
                PrimitiveGTO primitiveGTOxz(exponent,weight,pow9);
                PrimitiveGTO primitiveGTOyz(exponent,weight,pow10);


                contractedGTOxx.addPrimitive(primitiveGTOxx);
                contractedGTOyy.addPrimitive(primitiveGTOyy);
                contractedGTOzz.addPrimitive(primitiveGTOzz);
                contractedGTOxy.addPrimitive(primitiveGTOxy);
                contractedGTOxz.addPrimitive(primitiveGTOxz);
                contractedGTOyz.addPrimitive(primitiveGTOyz);
            }

            m_contractedGTOs.push_back(contractedGTOxx);
            m_contractedGTOs.push_back(contractedGTOyy);
            m_contractedGTOs.push_back(contractedGTOzz);
            m_contractedGTOs.push_back(contractedGTOxy);
            m_contractedGTOs.push_back(contractedGTOxz);
            m_contractedGTOs.push_back(contractedGTOyz);
        }

    }

}



rowvec BasisSet::corePosition() const
{
    return m_corePosition;
}

void BasisSet::setCorePosition(const rowvec &corePosition)
{
    m_corePosition = corePosition;
}

const ContractedGTO& BasisSet::getContracted(const int c) const
{
    return m_contractedGTOs.at(c);
}

int BasisSet::getNumContracted() const
{
    return m_contractedGTOs.size();
}

const int& BasisSet::getAngularMomentum() const
{
    return m_angularMomentum;
}

const int& BasisSet::coreCharge() const
{
    return m_coreCharge;
}

void BasisSet::setCoreCharge(const int &coreCharge)
{
    m_coreCharge = coreCharge;
}
const int& BasisSet::coreMass() const
{
    return m_coreMass;
}

void BasisSet::setCoreMass(const int &coreMass)
{
    m_coreMass = PROTONMASS * coreMass;
}





