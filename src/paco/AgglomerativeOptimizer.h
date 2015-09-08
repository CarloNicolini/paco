#ifndef _AGGLOMERATIVE_OPTIMIZER_H_
#define _AGGLOMERATIVE_OPTIMIZER_H_

#include "QualityOptimizer.h"

class AgglomerativeOptimizer : public QualityOptimizer
{
public:
    AgglomerativeOptimizer() {}
    AgglomerativeOptimizer(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb,const igraph_vector_t *weights=NULL);
    virtual ~AgglomerativeOptimizer();
    bool optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    void set_edges_order(const vector<int> &value);

protected:
    double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm);
    vector<int> edges_order;
};

#endif // AgglomerativeOptimizer_H
