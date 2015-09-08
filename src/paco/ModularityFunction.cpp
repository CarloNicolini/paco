#include "QualityFunction.h"
#include "ModularityFunction.h"
#include "igraph_utils.h"

ModularityFunction::ModularityFunction() {}

void ModularityFunction::eval(const igraph_t *g, igraph_vector_t *memb) const
{
    size_t n = igraph_vcount(g);
    size_t m = igraph_ecount(g);
    igraph_real_t density;
    igraph_density(g,&density,0);

    if (n != igraph_vector_size(memb) )
        throw std::runtime_error("Non consistent length of membership vector");

    // Initialize the vectors of edges and configuration model
    size_t nComms=(size_t)igraph_vector_max(memb)+1;

    igraph_vector_t observed;
    igraph_vector_init(&observed,nComms);

    // iterate all edges and check where the endpoints of the edges are
    // if they are in the same community
    igraph_integer_t v1,v2;
    for (igraph_integer_t edge_id=0; edge_id<m; edge_id++)
    {
        igraph_edge(g, edge_id, &v1, &v2);
        igraph_integer_t c1 = memb->stor_begin[v1]; // Community node "from" belongs
        igraph_integer_t c2 = memb->stor_begin[v2];  // Community node "to" belongs
        if (c1==c2)
            VECTOR(observed)[c1] += 1; // if in the same community, sum 1 because both endpoints of the vertex are in the same community.
    }
    //igraph_vector_print(&observed);
    // Sum the count of vertex pairs in every community
    quality = 0.0;
    for (igraph_integer_t c=0; c<nComms; c++)
    {
        size_t n_vertex_c = std::count(memb->stor_begin, memb->stor_end, c);
        double pairs_c = n_vertex_c *(n_vertex_c -1)/2;
        double density_c = VECTOR(observed)[c]/(pairs_c);
        quality += 2*pairs_c*KL(density_c,density);
    }

    igraph_vector_destroy(&observed);
}

void ModularityFunction::eval(const igraph_t *g, igraph_vector_t *memb, igraph_vector_t *weights) const
{
    size_t n = igraph_vcount(g);
    size_t m = igraph_vector_sum(weights);
    igraph_real_t density;
    igraph_density(g,&density,0);

    if (n != igraph_vector_size(memb) )
        throw std::runtime_error("Non consistent length of membership vector");

    // Initialize the vectors of edges and configuration model
    size_t nComms=(size_t)igraph_vector_max(memb)+1;

    igraph_vector_t observed;
    igraph_vector_init(&observed,nComms);

    // iterate all edges and check where the endpoints of the edges are
    // if they are in the same community
    igraph_integer_t v1,v2;
    for (igraph_integer_t edge_id=0; edge_id<m; edge_id++)
    {
        igraph_real_t w = weights->stor_begin[edge_id];
        igraph_edge(g, edge_id, &v1, &v2);
        igraph_integer_t c1 = memb->stor_begin[v1]; // Community node "from" belongs
        igraph_integer_t c2 = memb->stor_begin[v2];  // Community node "to" belongs
        if (c1==c2)
            VECTOR(observed)[c1] += w; // if in the same community, sum 1 because both endpoints of the vertex are in the same community.
    }
    //igraph_vector_print(&observed);
    // Sum the count of vertex pairs in every community
    quality = 0.0;
    for (igraph_integer_t c=0; c<nComms; c++)
    {
        size_t n_vertex_c = std::count(memb->stor_begin, memb->stor_end, c);
        double pairs_c = n_vertex_c *(n_vertex_c -1)/2;
        double density_c = VECTOR(observed)[c]/(pairs_c);
        quality += 2*pairs_c*KL(density_c,density);
    }

    igraph_vector_destroy(&observed);
}
