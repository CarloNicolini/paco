#ifndef _QUALITYOPTIMIZER_H
#define _QUALITYOPTIMIZER_H

#include "QualityFunction.h"
#include "PartitionHelper.h"
#include <set>

class QualityOptimizer
{
public:
    inline QualityOptimizer();
    inline QualityOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    inline virtual ~QualityOptimizer();
    virtual bool optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) = 0;
    const PartitionHelper* get_partition_helper() const;

protected:
    virtual double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm) = 0;
    PartitionHelper *par;
};

inline QualityOptimizer::QualityOptimizer() : par(NULL)
{
    par = new PartitionHelper();
}

inline QualityOptimizer::QualityOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights) : par(NULL)
{
    par = new PartitionHelper();
}

inline QualityOptimizer::~QualityOptimizer()
{
    if (par)
        delete par;
}

inline const PartitionHelper* QualityOptimizer::get_partition_helper() const
{
    return par;
}

#endif // _QUALITYOPTIMIZER_H

