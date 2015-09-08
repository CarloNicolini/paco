#ifndef _SIGNIFICANCEFUNCTION_H
#define _SIGNIFICANCEFUNCTION_H

#include "QualityFunction.h"
#include "KLDivergence.h"

class SignificanceFunction : public QualityFunction
{
public:
    SignificanceFunction();
    ~SignificanceFunction() {};

protected:
    void eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const;
    void eval(const PartitionHelper *par) const;
};

#endif // SignificanceFunction_H
