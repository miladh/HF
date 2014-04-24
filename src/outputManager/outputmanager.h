#ifndef OUTPUTMANAGER_H
#define OUTPUTMANAGER_H

#include <iostream>
#include <armadillo>
#include <boost/mpi.hpp>
#include <H5Cpp.h>

#include "../atom/atom.h"


using namespace arma;
using namespace std;
using namespace H5;

namespace hf{

class OutputManager
{
public:
    OutputManager(const int nAtoms, const string &outputFilePath);

    void saveAtoms(vector<Atom *> atoms);
    void saveEnergy(const double &energy, const mat &orbitalEnergies);
    void saveDipoleMoment(const double &dipoleMoment);
    void saveElectronDensity(const field<cube> &densityCubes);
    void closeOutput();

private:
    int m_rank;
    int m_nProcs;
    stringstream m_outputFileName;

    struct AtomAttributes {
        int type;
        char basisType[64];
        double x;
        double y;
        double z;
        int coreCharge;
        double corePartialCharge;
    };



    H5File *m_output;
    CompType *m_atomCompound;
    DataSet *m_dataset;
    AtomAttributes *m_atomAttributes;
};
}
#endif // OUTPUTMANAGER_H
