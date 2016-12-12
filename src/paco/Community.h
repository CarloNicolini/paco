/* This file is part of PACO-PArtitioning Clustering Optimization a program
* to find network partitions using modular solvers and quality functions.
*
*  Copyright (C) 2016 Carlo Nicolini <carlo.nicolini@iit.it>
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
*
* If you use PACO for you publication please cite:
*
* "Modular structure of brain functional networks: breaking the resolution limit by Surprise"
* C. Nicolini and A. Bifone, Scientific Reports
* doi:10.1038/srep19250
*
* "Community detection in weighted brain connectivity networks beyond the resolution limit", 
* C.Nicolini. C.Bordier, A.Bifone, arxiv 1609.04316
* https://arxiv.org/abs/1609.04316
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

enum OptimizerType
{
    MethodAgglomerative = 0,
    MethodRandom = 1,
    MethodAnneal = 2,
    MethodInfomap = 3
};

enum QualityType
{
    QualitySurprise = 0,
    QualitySignificance = 1,
    QualityAsymptoticSurprise = 2,
    QualityInfoMap = 3,
#ifdef EXPERIMENTAL
    QualityModularity = 4,
    QualityAsymptoticModularity = 5,
    QualityWonder = 6,
    QualityDegreeCorrectedSurprise = 7
#endif
};

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
    vector<int> get_membership_vector() const; // for python version of PACO
    size_t get_membership(size_t i) const;
    void reindex_membership();
    void order_membership();
    void print_membership();
    void save_membership(const char *filename, const igraph_vector_t *m=NULL);

    // Methods needed by agglomerative optimizer
    void sort_edges();
    vector<int> get_sorted_edges_indices();
    double optimize(QualityType qual, OptimizerType opt, int nrep=1);

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
