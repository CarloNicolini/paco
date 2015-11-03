/* This file is part of PACO-PArtitioning Clustering Optimization a program
* to find network partitions using modular solvers and quality functions.
*
*  Copyright (C) 2015 Carlo Nicolini <carlo.nicolini@iit.it>
*
*  PACO is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 3 of the License, or (at your option) any later version.
*
*  Alternatively, you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  PACO is distributed in the hope that it will be useful, but WITHOUT ANY
*  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
*  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General Public
*  License and a copy of the GNU General Public License along with
*  PACO. If not, see <http://www.gnu.org/licenses/>.
*/

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
        {
            return post-pre;
        }
    }
    else
        return 0;
}

void AgglomerativeOptimizer::set_edges_order(const vector<int> &value)
{
    edges_order = value;
}

double AgglomerativeOptimizer::optimize(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb,const igraph_vector_t *weights)
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
    //par->reindex(memb);
    return fun(g,memb,weights);
#ifdef _DEBUG
    par->print();
    printf(ANSI_COLOR_RED "AGGLOMERATIVE Final Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
}
