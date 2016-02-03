#include "QualityFunction.h"
#include "CorrelationModularityFunction.h"
#include "KLDivergence.h"
#include "igraph_utils.h"

CorrelationModularityFunction::CorrelationModularityFunction()
{

}

void CorrelationModularityFunction::eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights) const
{
    PartitionHelper *par2 = new PartitionHelper;
    par2->init(g,memb,weights);
    this->eval(par2);
    delete par2;
}

void CorrelationModularityFunction::eval(const PartitionHelper *par) const
{
    quality = 0;
    igraph_real_t m = par->get_graph_total_weight();

    if (m > 0)
    {
        for (CommMapCIter iter=par->get_communities().begin(); iter!=par->get_communities().end(); ++iter)
        {
            size_t c = iter->first;
            igraph_real_t a = par->get_incomm_weight().at(c)/m;
            igraph_real_t e = par->get_incomm_deg().at(c)/(2*m);
            quality += KL(a,e*e);
        }
    }
}
