#ifndef SUMFUNCTION_H
#define SUMFUNCTION_H

#include "QualityFunctionImpl.h"
#include "UnionFunction.h"
#include <iostream>

class SumFunction : public UnionFunction
{
  public:

    SumFunction(QualityFunctionImpl* f1, QualityFunctionImpl* f2 )
      : UnionFunction(f1,f2){}

    ~SumFunction(){}

    void eval(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const
    {
      std::cout << "(";
      quality = (*f1)(g,m,weights);
      std::cout << "+";
      quality += (*f2)(g,m,weights);
      std::cout << ")";
    }

    QualityFunctionImpl* clone() const
    {
        return new SumFunction(f1->clone(),f2->clone());
    }
};

#endif // SUMFUNCTION_H

