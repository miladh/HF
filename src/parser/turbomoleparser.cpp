#include "turbomoleparser.h"

using namespace hf;

TurbomoleParser::TurbomoleParser(string filename)
{
    loadfile(filename);
}

void TurbomoleParser::loadfile(string filename)
{

    string line, stringToSearch;
    fstream file;
    file.open(filename);


    //-----------------------------------------------------------------------------------
    //Read turbomole file
    if (file.is_open()){
        while (getline(file, line)){
            stringToSearch += line +"\n";
        }
        file.close();
    } else {
        cout << "Error: Could not open file "<< filename << endl;
        exit(EXIT_FAILURE);
    }

    //-----------------------------------------------------------------------------------
    regex typeRegex("\\s*([a-zA-Z]+)\\s*(([0-9])-([0-9]+)([G]))\\s*");
    sregex_iterator type(stringToSearch.begin(), stringToSearch.end(), typeRegex);
    sregex_iterator endType;

    for(; type!=endType; type++){
        //Search atom type and basis
        bool skip = false;
        skip |= regex_match(stringToSearch, regex("#.*"));
        skip |= regex_match(stringToSearch, regex("$basis.*"));
        skip |= regex_match(stringToSearch, regex("$end.*"));
        skip |= regex_match(stringToSearch, regex("\\*.*"));
        if(skip) {
            continue;
        }
        smatch what;
        while(regex_search(stringToSearch, what, typeRegex)) {
            AtomMeta atomMeta = AtomMeta::getData(string(what[1]));
            setAtomType(atomMeta.type);
            setAtomCharge(atomMeta.charge);
            setAtomMass(atomMeta.mass);
            break;
        }

    }


    //-----------------------------------------------------------------------------------
    //Search for groups of contracted GTOs
    regex searchCGTO("([0-9]\\s+[spd])((\\s+-?[0-9]+\\.[0-9]+)+)");
    sregex_iterator CGTOs(stringToSearch.begin(), stringToSearch.end(), searchCGTO);
    sregex_iterator endCGTOs;

    for(; CGTOs!=endCGTOs; CGTOs++){
        string subString = CGTOs->str(2).c_str();

        //Search for groups of (exponent and coefficient) for each primitive
        regex searchPGTO("(-?[0-9]+\\.[0-9]+)\\s+(-?[0-9]+\\.[0-9]+)\\s*");
        sregex_iterator PGTOs(subString.begin(), subString.end(), searchPGTO);
        sregex_iterator endPGTOs;

        // Check type of orbital (s, p, d)
        string orbitalType = CGTOs->str(1).c_str();
        regex searchs("[0-9]\\s+s");
        regex searchp("[0-9]\\s+p");
        regex searchd("[0-9]\\s+d");



        // s-type:
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

        // p-type:
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

        // d-type:
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
int TurbomoleParser::angularMomentum() const
{
    return m_angularMomentum;
}

void TurbomoleParser::setAngularMomentum(int angularMomentum)
{
    m_angularMomentum = angularMomentum;
}

double TurbomoleParser::atomMass() const
{
    return m_atomMass;
}

void TurbomoleParser::setAtomMass(double atomMass)
{
    m_atomMass = atomMass;
}

int TurbomoleParser::atomCharge() const
{
    return m_atomCharge;
}

void TurbomoleParser::setAtomCharge(int atomCharge)
{
    m_atomCharge = atomCharge;
}


vector<ContractedGTO> TurbomoleParser::contractedGTOs() const
{
    return m_contractedGTOs;
}

int TurbomoleParser::atomType() const
{
    return m_atomType;
}

void TurbomoleParser::setAtomType(int atomType)
{
    m_atomType = atomType;
}



