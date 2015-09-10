#ifndef SURPRISEFUNCTIONPIMPL_H
#define SURPRISEFUNCTIONPIMPL_H

#include "QualityFunctionImpl.h"

class SurpriseFunction2 : public QualityFunctionImpl
{
public:
    SurpriseFunction2();
    ~SurpriseFunction2(){}

protected:
    void eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL) const;
    void eval(const PartitionHelper *par) const;
    QualityFunctionImpl* clone() const;
};

#endif // SURPRISEFUNCTIONPIMPL_H
