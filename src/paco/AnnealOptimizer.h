#ifndef _ANNEALOPTIMIZER_H
#define _ANNEALOPTIMIZER_H

#include "QualityOptimizer.h"

class AnnealParameters
{
public:
    AnnealParameters(size_t nIter=1E3, size_t maxHits=1E2, double mintemp=1E-1, double tol=1E-5, double temp=1E6, double tempscale=0.99, double minfval=std::numeric_limits<double>::max()) :
        nIterations(nIter),
        nHits(maxHits),
        min_temp(mintemp),
        tolerance(tol),
        temperature(temp),
        temp_scale(tempscale),
        minfval(minfval)
    {}
    size_t nIterations;
    size_t nHits;
    double min_temp;
    double tolerance;
    double temperature;
    double temp_scale;
    double minfval;
};

class AnnealOptimizer : public QualityOptimizer
{
public:
    AnnealOptimizer() {}
    AnnealOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    virtual ~AnnealOptimizer();
    void set_parameters(const AnnealParameters &_par)
    {
        this->param = _par;
    }
    bool optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);

protected:
    double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm);

    AnnealParameters param;
};
#endif // _ANNEALOPTIMIZER_H
