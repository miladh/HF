#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <H5Cpp.h>

#include "../atom/atom.h"


using namespace arma;
using namespace std;

namespace hf{

class OutputManager
{
public:
    OutputManager();

    void saveEnergy(const double &energy);
    void registerAtoms(vector<Atom *> atoms);

private:
    int m_rank;
    int m_nProcs;
    stringstream m_outputFileName;
    H5::H5File *m_output;

    struct AtomProperties {
        int type;
        char basisType[64];
        double x;
        double y;
        double z;
        double corePartialCharge;
        int coreCharge;
    };


};
}
#endif // OUTPUTMANAGER_H
