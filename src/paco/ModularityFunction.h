#ifndef _MODULARITYFUNCTION_H
#define _MODULARITYFUNCTION_H

#include "QualityFunction.h"
#include "KLDivergence.h"

class ModularityFunction : public QualityFunction
{
public:
    ModularityFunction();
    ~ModularityFunction() {};

protected:
    void eval(const igraph_t *g, igraph_vector_t *memb) const;
    void eval(const igraph_t *g, igraph_vector_t *memb, igraph_vector_t *weights) const;
};

#endif // _MODULARITYFUNCTION_H
