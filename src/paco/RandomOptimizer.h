#ifndef _RANDOPTIMIZER_H_
#define _RANDOPTIMIZER_H_

#include "QualityOptimizer.h"

class RandomOptimizer : public QualityOptimizer
{
public:
    RandomOptimizer() {}
    RandomOptimizer(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    virtual ~RandomOptimizer();
    bool optimize(const igraph_t *g, const QualityFunction &fun,const  igraph_vector_t *memb, const igraph_vector_t *weights=NULL);

protected:
    double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm);

};

#endif // RANDOMOPTIMIZER_H
