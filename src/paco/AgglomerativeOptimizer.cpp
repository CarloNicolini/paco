#include "AgglomerativeOptimizer.h"
#include <iostream>

AgglomerativeOptimizer::AgglomerativeOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights) : QualityOptimizer(g, fun, memb)
{
    if (edges_order.empty())
    {
        for (size_t i=0; i<par->get_num_edges(); ++i)
            edges_order.push_back(i);
    }
    this->optimize(g,fun,memb,weights);
}

AgglomerativeOptimizer::~AgglomerativeOptimizer()
{
}

double AgglomerativeOptimizer::diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm)
{
    // try to join the vertices
    double pre = fun(par);
    int orig_comm = memb->stor_begin[vert]; // save old original community of vert
    bool vertex_moved = par->move_vertex(g, memb,vert,dest_comm);
    if (vertex_moved)
    {
        double post = fun(par);
        if (post<pre)
        {
            par->move_vertex(g, memb,vert,orig_comm); // restore previous status
            return post-pre;
        }
        else
            return post-pre;
    }
    else
        return 0;
}

void AgglomerativeOptimizer::set_edges_order(const vector<int> &value)
{
    edges_order = value;
}

bool AgglomerativeOptimizer::optimize(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb,const igraph_vector_t *weights)
{
    par->init(g,memb,weights);
    if (edges_order.empty())
    {
        for (size_t i=0; i<par->get_num_edges(); ++i)
            edges_order.push_back(i);
    }

#ifdef _DEBUG
    printf(ANSI_COLOR_RED "AGGLOMERATIVE Initial Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
    size_t m = edges_order.size();

    for (size_t i=0; i<m; ++i)
    {
        int e = edges_order.at(i); // edge to consider
        int vert1, vert2;
        igraph_edge(g,e,&vert1,&vert2);
        if ( rand()%2 ) // Randomly choose to aggregate vert1-->comm2 or vert2-->comm1
        {
            size_t dest_comm = memb->stor_begin[vert2];
            diff_move(g,fun,memb,vert1,dest_comm);
        }
        else
        {
            size_t dest_comm = memb->stor_begin[vert1];
            diff_move(g,fun,memb,vert2,dest_comm);
        }
    }
    par->reindex(memb);
#ifdef _DEBUG
    par->print();
    printf(ANSI_COLOR_RED "AGGLOMERATIVE Final Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
}
