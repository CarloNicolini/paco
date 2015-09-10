#ifndef QUALITYFUNCTIONIMPL_H
#define QUALITYFUNCTIONIMPL_H

#include <iostream>

#include <igraph.h>
#include "igraph_utils.h"

#include "PartitionHelper.h"

using std::cout;
using std::cerr;
using std::endl;

class QualityFunctionImpl
{
public:

    virtual ~QualityFunctionImpl() {}

    double &operator()()
    {
        return quality;
    }

    double &operator()(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const
    {
        eval(g,m,weights);
        return quality;
    }

    double &operator()(const PartitionHelper *par) const
    {
        eval(par);
        return quality;
    }

    virtual QualityFunctionImpl* clone() const = 0;

protected:
    mutable double quality; // mutable keyword allows to modify quality even in const methods.
    virtual void eval(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const = 0;
    virtual void eval(const PartitionHelper *par) const = 0;
    PartitionHelper *par;

};

#endif // QUALITYFUNCTIONIMPL_H

