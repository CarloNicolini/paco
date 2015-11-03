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

#include "QualityFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "igraph_utils.h"
#include "KLDivergence.h"

/**
 * @brief AsymptoticSurpriseFunction::AsymptoticSurpriseFunction
 */
AsymptoticSurpriseFunction::AsymptoticSurpriseFunction() {}

/**
 * @brief AsymptoticSurpriseFunction::eval
 * @param g
 * @param memb
 * @param weights
 */
void AsymptoticSurpriseFunction::eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights) const
{
    /*
    igraph_real_t n = static_cast<igraph_real_t>(igraph_vcount(g));
    igraph_real_t num_edges = static_cast<igraph_real_t>(igraph_ecount(g));
    igraph_real_t m = igraph_real_t(igraph_ecount(g));
    igraph_real_t p = n*(n-1)/2;

    if (weights)
        m = igraph_vector_sum(weights);

    if (n != igraph_vector_size(memb) )
        throw std::runtime_error("Non consistent length of membership vector");

    // Sum of intracluster edge weights
    igraph_real_t mi=0;
    // Sum of intracluster pairs
    igraph_real_t pi=0;

    // Initialize the vectors of edges and configuration model
    size_t nComms=(size_t)igraph_vector_max(memb)+1; // XXX to fix in a future...

    igraph_integer_t from;
    igraph_integer_t to;

    // iterate all edges and check where the endpoints of the edges are
    // if they are in the same community
    for (igraph_integer_t edge_id=0; edge_id<num_edges; edge_id++)
    {
        igraph_edge(g, (igraph_integer_t) edge_id, &from, &to);
        igraph_real_t w = 1;
        if (weights)
            w = weights->stor_begin[edge_id];
        igraph_integer_t comm_from=*(memb->stor_begin+from); // Community node "from" belongs
        igraph_integer_t comm_to=*(memb->stor_begin+to);  // Community node "to" belongs
        int delta = static_cast<igraph_real_t>(int(comm_from==comm_to));
        mi += delta*w;
    }

    // Sum the count of vertex pairs in every community
    for (igraph_integer_t c=0; c<nComms; c++)
    {
        igraph_real_t vertices_count = static_cast<igraph_real_t>(std::count(memb->stor_begin, memb->stor_end, c));
        pi += vertices_count*(vertices_count-1.0)/2.0;
    }

    double observed = mi/m;
    double expected = pi/p;

    quality = m*KL(observed,expected);
    */
    PartitionHelper *par = new PartitionHelper;
    par->init(g,memb,weights);
    this->eval(par);
    delete par;
    //printf("--> AS=%f -- mi=%f pi=%f m=%f p=%f\n",quality, mi,pi,m,p);
}

/**
 * @brief AsymptoticSurpriseFunction::eval
 * @param par
 */
void AsymptoticSurpriseFunction::eval(const PartitionHelper *par) const
{
    /*
    double p = par->get_graph_total_pairs();
    double pi = par->get_total_incomm_pairs();
    double m = par->get_graph_total_weight();
    double mi = par->get_total_incomm_weight();

    double observed = mi/m;
    double expected = pi/p;
    quality = m*KL(observed,expected);
    */
    //printf("AS=%f -- mi=%f pi=%f m=%f p=%f\n",quality, mi,pi,m,p);
    quality = 0;
    igraph_real_t m = par->get_graph_total_weight();
    igraph_real_t p = par->get_graph_total_pairs();

    igraph_real_t mi = 0;
    igraph_real_t pi = 0;

    if (m > 0)
    {
        for (CommMapCIter iter=par->get_communities().begin(); iter!=par->get_communities().end(); ++iter)
        {
            size_t c = iter->first;
            igraph_real_t mc = par->get_incomm_weight().at(c);
            igraph_real_t nc = par->get_incomm_nvert().at(c);
            igraph_real_t pc = nc*(nc-1.0)/2.0;
            mi += mc;
            pi += pc;
        }
    }
    quality = m*KL(mi/m,pi/p);
}
