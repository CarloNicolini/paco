#ifndef AsymptoticSurpriseFunction_H
#define AsymptoticSurpriseFunction_H

#include "QualityFunction.h"
#include "Surprise.h"

class AsymptoticSurpriseFunction : public QualityFunction
{
public:
    AsymptoticSurpriseFunction();
    ~AsymptoticSurpriseFunction() {};

protected:
    void eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const;
    void eval(const PartitionHelper *par) const;
};

#endif // AsymptoticSurpriseFunction_H
