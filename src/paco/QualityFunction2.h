#ifndef QUALITYFUNCTION2_H
#define QUALITYFUNCTION2_H

#include "QualityFunctionImpl.h"
#include "SumFunction.h"
#include "UnionFunction.h"

class QualityFunction2
{
  public:

    QualityFunction2( QualityFunctionImpl* impl )
      : impl(impl)
    { }

    ~QualityFunction2()
    {
        delete impl;
    }

    double &operator()()
    {
      return (*impl)();
    }

    double &operator()(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const
    {
      return (*impl)(g,m,weights);
    }

    // This method can be friend
    QualityFunction2 operator+(const QualityFunction2& f2) const
    {
      return QualityFunction2(new SumFunction(impl->clone(), f2.impl->clone()));
    }

  private:
    QualityFunctionImpl* impl;
};

#endif // QUALITYFUNCTION2_H

