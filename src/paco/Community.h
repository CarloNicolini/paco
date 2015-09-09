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

#ifndef _COMMUNITY_H_
#define _COMMUNITY_H_

#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>

#include "Graph.h"

typedef std::map<igraph_integer_t, igraph_vs_t> CommunityVertexMap;
typedef std::map<igraph_integer_t, igraph_vs_t>::iterator CommunityVertexMapIterator;

class CommunityStructure
{
public:
    CommunityStructure(const GraphC *pgraph);
    ~CommunityStructure();
    void init(const GraphC *pgraph);
    void set_random_seed(int seed=-1); // XXX da spostare...

    // Methods that work on membership
    void read_membership_from_file(const std::string &filename);
    const igraph_vector_t* get_membership() const;
    size_t get_membership(size_t i) const;
    void reindex_membership();
    void print_membership();
    void save_membership(const char *filename, const igraph_vector_t *m=NULL);

    // Methods needed by agglomerative optimizer
    void sort_edges();
    vector<int> get_sorted_edges_indices();

protected:
    void compute_pairwise_similarities();
    void compute_edges_similarities();

private:
    const GraphC* pgraph; // internal pointer to Graph proxy
    igraph_integer_t nVertices;
    igraph_integer_t nEdges;
    igraph_integer_t nPairs;

    igraph_vector_t membership;

    // Random number generator
    igraph_rng_t rng;

    // Needed by agglomerative optimizer
    igraph_vector_t edges_sim;
    std::vector < std::pair<int, igraph_real_t> > sorted_edges;
};

#endif
