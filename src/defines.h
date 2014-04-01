#ifndef DEFINES_H
#define DEFINES_H

#include <boost/regex.hpp>
#include <boost/config.hpp>

#define HFSOLVERTOLERANCE 1.0E-8
#define PROTONMASS 1836
//#define h 1.0E-5


#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

//#define USE_MPI
struct mpiTask{
    int nFunctionCalls;
    bool isAvailable;
    int p;
};




#endif // DEFINES_H
