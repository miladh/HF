#ifndef ANALYSER_H
#define ANALYSER_H

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;

#include<system/system.h>
#include<hfSolver/hfsolver.h>


class Analyser
{
public:
    Analyser(System *system, int nGridPoints);
    void calculatePES();
    void calculateForces();
private:
    System *m_system;
    vec bondLength, m_energy, m_Force;
    HFsolver* solver;

    void writeToFile(const vec R, int n);
};

#endif // ANALYSER_H
