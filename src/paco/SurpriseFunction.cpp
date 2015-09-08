#include "QualityFunction.h"
#include "SurpriseFunction.h"
#include "igraph_utils.h"

SurpriseFunction::SurpriseFunction() {}

void SurpriseFunction::eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights) const
{
    size_t n = igraph_vcount(g);
    size_t m = igraph_ecount(g);
    size_t p = n*(n-1)/2;

    if (n != igraph_vector_size(memb) )
        throw std::runtime_error("Non consistent length of membership vector");

    // Sum of intracluster edge weights
    size_t mzeta=0;
    // Sum of intracluster pairs
    size_t pzeta=0;

    // Initialize the vectors of edges and configuration model
    size_t nComms=(size_t)igraph_vector_max(memb)+1; // XXX to fix in a future...

    igraph_integer_t from;
    igraph_integer_t to;

    // iterate all edges and check where the endpoints of the edges are
    // if they are in the same community
    for (igraph_integer_t edge_id=0; edge_id<m; edge_id++)
    {
        igraph_edge(g, (igraph_integer_t) edge_id, &from, &to);
        igraph_integer_t comm_from=*(memb->stor_begin+from); // Community node "from" belongs
        igraph_integer_t comm_to=*(memb->stor_begin+to);  // Community node "to" belongs
        mzeta += size_t(comm_from==comm_to);
    }

    // Sum the count of vertex pairs in every community
    for (igraph_integer_t c=0; c<nComms; c++)
    {
        size_t vertices_count = std::count(memb->stor_begin, memb->stor_end, c);
        pzeta += vertices_count*(vertices_count-1)/2;
    }
    quality = computeSurprise(p,pzeta,m,mzeta);
}

void SurpriseFunction::eval(const PartitionHelper *par) const
{
    double p = par->get_graph_total_pairs();
    double pi = par->get_total_incomm_pairs();
    double m = par->get_graph_total_weight();
    double mi = par->get_total_incomm_weight();

    quality = computeSurprise(p,pi,m,mi);
}
