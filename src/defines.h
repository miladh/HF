#ifndef DEFINES_H
#define DEFINES_H

#include <boost/regex.hpp>
#include <boost/config.hpp>

#define HFSOLVERTOLERANCE 1.0E-8
#define PROTONMASS 1836
#define h 1.0E-5

#define USE_MPI 1
struct mpiTask{
    int nFunctionCalls;
    bool isAvailable;
    int p;
};

#endif // DEFINES_H
