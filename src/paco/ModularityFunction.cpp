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
