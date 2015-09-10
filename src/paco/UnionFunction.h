#ifndef UNIONFUNCTION_H
#define UNIONFUNCTION_H

#include "QualityFunctionImpl.h"

class UnionFunction : public QualityFunctionImpl
{
  public:

    UnionFunction( QualityFunctionImpl* f1, QualityFunctionImpl* f2 )
      : f1(f1), f2(f2)
    { }

    ~UnionFunction()
    { delete f1; delete f2; }

  //protected:
    QualityFunctionImpl* f1;
    QualityFunctionImpl* f2;
};
#endif // UNIONFUNCTION_H

