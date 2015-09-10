#ifndef SUMFUNCTION_H
#define SUMFUNCTION_H

#include <iostream>
#include "QualityFunctionImpl.h"
#include "UnionFunction.h"

class SumFunction : public UnionFunction
{
  public:

    SumFunction(QualityFunctionImpl* f1, QualityFunctionImpl* f2 )
      : UnionFunction(f1,f2){}

    ~SumFunction(){}

    void eval(const igraph_t *g, const igraph_vector_t *m, const igraph_vector_t *weights=NULL) const
    {
      std::cout << "( ";
      quality = (*f1)(g,m,weights);
      std::cout << quality << " + ";
      quality += (*f2)(g,m,weights);
      std::cout << " )";
    }

    void eval(const PartitionHelper *par) const
    {
        quality = f1->operator ()(par) + f2->operator ()(par);
    }

    QualityFunctionImpl* clone() const
    {
        return new SumFunction(f1->clone(),f2->clone());
    }
};

#endif // SUMFUNCTION_H

