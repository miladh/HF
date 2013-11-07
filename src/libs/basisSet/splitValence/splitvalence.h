#ifndef SPLITVALENCE_H
#define SPLITVALENCE_H

#include<basisSet/basisset.h>

class SplitValence : public BasisSet
{
public:
    SplitValence();

protected:
    void normalizeGaussian();

};

#endif // SPLITVALENCE_H
