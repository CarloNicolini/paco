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

AsymptoticSurpriseFunction::AsymptoticSurpriseFunction() {}

void AsymptoticSurpriseFunction::eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights) const
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
    quality = m*KL(double(mzeta)/double(m),double(pzeta)/double(p));
}

void AsymptoticSurpriseFunction::eval(const PartitionHelper *par) const
{
    double p = par->get_graph_total_pairs();
    double pi = par->get_total_incomm_pairs();
    double m = par->get_graph_total_weight();
    double mi = par->get_total_incomm_weight();

    quality = m*KL(mi/m,pi/p);
}
