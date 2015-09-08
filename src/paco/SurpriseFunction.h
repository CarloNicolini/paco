#ifndef SURPRISEFUNCTION_H
#define SURPRISEFUNCTION_H

#include "QualityFunction.h"
#include "Surprise.h"

class SurpriseFunction : public QualityFunction
{
public:
    SurpriseFunction();
    ~SurpriseFunction() {};

protected:
    void eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const;
    void eval(const PartitionHelper *par) const;
};

#endif // SURPRISEFUNCTION_H
