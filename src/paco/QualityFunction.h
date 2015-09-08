#ifndef QUALITYFUNCTION
#define QUALITYFUNCTION

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <memory>

#include <igraph.h>
#include "igraph_utils.h"

#include "PartitionHelper.h"

using std::cout;
using std::cerr;
using std::endl;

class QualityFunction
{
protected:
    mutable double quality; // mutable keyword allows to modify quality even in const methods.
    // Takes const pointer to igraph_t so that graph is unmodifiable, to implement in children classes.
    virtual void eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const = 0;
    virtual void eval(const PartitionHelper *par) const = 0;

public:
    virtual ~QualityFunction() {}

    double &operator()()
    {
        return quality;
    }

    double &operator()(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const
    {
        eval(g,memb,weights);
        return quality;
    }

    double &operator()(const PartitionHelper *par) const
    {
        eval(par);
        return quality;
    }
};


#endif // QUALITYFUNCTION
