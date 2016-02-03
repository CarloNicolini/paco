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

/**
 * @brief AgglomerativeOptimizer::AgglomerativeOptimizer
 * @param g
 * @param fun
 * @param memb
 * @param weights
 */
AgglomerativeOptimizer::AgglomerativeOptimizer(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights) : QualityOptimizer(g, fun, memb)
{
    if (edges_order.empty())
    {
        for (size_t i=0; i<par->get_num_edges(); ++i)
            edges_order.push_back(i);
    }
    this->optimize(g,fun,memb,weights);
}

/**
 * @brief AgglomerativeOptimizer::~AgglomerativeOptimizer
 */
AgglomerativeOptimizer::~AgglomerativeOptimizer()
{
}

/**
 * @brief AgglomerativeOptimizer::diff_move
 * @param g
 * @param fun
 * @param memb
 * @param vert
 * @param dest_comm
 * @return
 */
#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include <typeinfo>

double AgglomerativeOptimizer::diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm, const igraph_vector_t *weights)
{
    int orig_comm = memb->stor_begin[vert]; // save old original community of vert

    // Control the type of quality function with RTTI
    // Heuristic to skip edges in the same community. Surprise and AsymptoticSurprise do not change if the edge is already intracluster
    if (orig_comm == (int)dest_comm && (typeid(fun)==typeid(SurpriseFunction) || typeid(fun)==typeid(AsymptoticSurpriseFunction)) )
        return 0.0;

    // try to join the vertices
    double pre = fun(par);
    bool vertex_moved = par->move_vertex(g, memb,vert,dest_comm,weights);


    if (vertex_moved)
    {
        double post = fun(par);
        if (post<pre)
        {
            par->move_vertex(g, memb,vert,orig_comm,weights); // restore previous status
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

/**
 * @brief AgglomerativeOptimizer::set_edges_order
 * @param value
 */
void AgglomerativeOptimizer::set_edges_order(const vector<int> &value)
{
    edges_order = value;
}

/**
 * @brief AgglomerativeOptimizer::optimize
 * @param g
 * @param fun
 * @param memb
 * @param weights
 * @return
 */
double AgglomerativeOptimizer::optimize(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb, const igraph_vector_t *weights)
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
        double deltaS=0;
        //printf(ANSI_COLOR_RED "Evaluating edge %d-%d\n",vert1,vert2);

        if ( rand()%2 ) // Randomly choose to aggregate vert1-->comm2 or vert2-->comm1
        {
            size_t dest_comm = memb->stor_begin[vert2];
            deltaS=diff_move(g,fun,memb,vert1,dest_comm,weights);
        }
        else
        {
            size_t dest_comm = memb->stor_begin[vert1];
            deltaS=diff_move(g,fun,memb,vert2,dest_comm,weights);
        }
        /*
        if (deltaS>0)
        {
            printf(ANSI_COLOR_RED "Edge %d-%d merged\n",vert1,vert2);
            par->print();
            printf(ANSI_COLOR_GREEN "Surprise increment %g\n",deltaS);
            printf(ANSI_COLOR_GREEN "==============\n");
        }
        if (deltaS<=0)
        {
            printf(ANSI_COLOR_RED "Edge %d-%d NOT merged\n",vert1,vert2);
            printf(ANSI_COLOR_GREEN "Surprise increment %g\n",deltaS);
            printf(ANSI_COLOR_GREEN "==============\n");
        }
        */
    }
    //par->reindex(memb);
    //par->print();
    return fun(par);

#ifdef _DEBUG
    par->print();
    printf(ANSI_COLOR_RED "AGGLOMERATIVE Final Qual=%g\n" ANSI_COLOR_RESET,fun(par));
#endif
}
