#include "RandomOptimizer.h"
#include <iostream>

RandomOptimizer::RandomOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights) : QualityOptimizer(g, fun, memb)
{
    this->optimize(g,fun,memb,weights);
}

RandomOptimizer::~RandomOptimizer()
{
}

double RandomOptimizer::diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm)
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

bool RandomOptimizer::optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights)
{
    par->init(g,memb,weights);
#ifdef _DEBUG
    printf(ANSI_COLOR_RED "RANDOM Initial Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
    for (int i=0; i<igraph_ecount(g); i++)
    {
        int e = rand()%igraph_ecount(g);
        int vert1;
        int vert2;
        igraph_edge(g,e,&vert1,&vert2);

        size_t dest_comm = memb->stor_begin[vert2];

        diff_move(g,fun,memb,vert1,dest_comm);
    }
    par->reindex(memb);
    #ifdef _DEBUG
    par->print();
    printf(ANSI_COLOR_RED "RANDOM Final Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
}
