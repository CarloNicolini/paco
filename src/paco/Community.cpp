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

#include <igraph.h>
#include "Community.h"

/**
 * @brief CommunityStructure::CommunityStructure
 * @param G
 */
CommunityStructure::CommunityStructure(const GraphC *G)
{
    this->init(G);
}

/**
 * @brief CommunityStructure::~CommunityStructure
 */
CommunityStructure::~CommunityStructure()
{
    igraph_vector_destroy(&this->membership);
    igraph_vector_destroy(&edges_sim);
    igraph_rng_destroy(&this->rng);
}

/**
 * @brief CommunityStructure::read_membership_from_file
 * @param filename
 */
void CommunityStructure::read_membership_from_file(const string &filename)
{
    std::ifstream membership_file(filename.c_str());
    if (!membership_file.good())
        throw std::runtime_error(std::string("File not found"));

    string line;
    vector<igraph_real_t> vmemb;
    while ( getline(membership_file,line) )
    {
        if (line.size() == 0)
            continue;
        std::stringstream str(line);
        igraph_real_t memb;
        str >> memb;
        vmemb.push_back(memb);
    }
    if (vmemb.size() != this->nVertices)
        throw std::runtime_error("Error loading membership file. Vertex number not consistent with current graph.");

    for (size_t i=0; i<this->nVertices; ++i)
        VECTOR(membership)[i] = vmemb[i];
}

void CommunityStructure::save_membership(const char *filename, const igraph_vector_t *m)
{
    if (m==NULL)
    {
        save_membership(filename,&membership);
    }
    else
    {
        std::ofstream outputmembership;
        outputmembership.open(filename);
        if (!outputmembership.good())
            throw std::ios_base::failure("Error, file" + std::string(filename)+ " exists");
        for (size_t i=0; i<igraph_vector_size(m); i++)
            outputmembership << m->stor_begin[i] << endl;
        outputmembership.close();
    }
}

/**
 * @brief CommunityStructure::init
 * @param G
 */
void CommunityStructure::init(const GraphC *G)
{
    // Copy the graph instance pointer
    this->pgraph = G;

    this->nVertices = G->number_of_nodes();
    this->nEdges = G->number_of_edges();
    this->nPairs = nVertices*(nVertices-1)/2;

    // Init the membership vector
    IGRAPH_TRY(igraph_vector_init(&membership,nVertices));
    for (size_t i=0; i<nVertices; ++i)
        VECTOR(membership)[i]=i;

    // Initialize the vector containing the edges similarities
    IGRAPH_TRY(igraph_vector_init(&edges_sim,0));

    // Initialize the random number generator
    IGRAPH_TRY(igraph_rng_init(&this->rng,&igraph_rngtype_mt19937));
}

void CommunityStructure::print_membership()
{
    igraph_vector_print(&membership);
}

/**
 * @brief CommunityStructure::set_random_seed, if negative, then seed based on current time is used (default seed=-1)
 * @param seed
 */
void CommunityStructure::set_random_seed(int seed)
{
    if (seed<0)
    {
#ifdef WIN32
        QueryPerformanceCounter(&endCount);
        std::srand(startCount.QuadPart);
#endif
#if defined(__linux__) || defined(__APPLE__) 
        struct timeval start;
        gettimeofday(&start, NULL);
        std::srand(start.tv_usec);
#endif
        IGRAPH_TRY(igraph_rng_seed(&this->rng,start.tv_usec));
    }
    else
    {
        std::srand(size_t(seed));
        IGRAPH_TRY(igraph_rng_seed(&this->rng,seed));
    }
}

/**
 * @brief CommunityStructure::reindex_membership
 */
void CommunityStructure::reindex_membership()
{
    IGRAPH_TRY(igraph_reindex_membership(&this->membership,NULL));
    int minC = igraph_vector_min(&this->membership);
    for (size_t i=0; i<nVertices; ++i)
        VECTOR(membership)[i] -= minC;
}

/**
 * @brief CommunityStructure::compute_edges_similarities
 */
void CommunityStructure::compute_edges_similarities()
{
    // Query the similarities for all edges
    igraph_similarity_jaccard_es(this->pgraph->get_igraph(), &edges_sim, igraph_ess_all(IGRAPH_EDGEORDER_ID), IGRAPH_ALL, false);
}

/**
 * @brief CommunityStructure::sort_edges
 */
void CommunityStructure::sort_edges()
{
    this->compute_edges_similarities(); // initialize the edges_sim igraph_vector_t
    sorted_edges.resize(nEdges);

    for (size_t i=0; i<this->nEdges; ++i)
    {
        sorted_edges.at(i).first = i;
        sorted_edges.at(i).second = edges_sim.stor_begin[i];
    }

    // Edges are sorted by jaccard index. Edges endpoints can be found by IGRAPH_TO and IGRAPH_FROM macros
    std::sort(sorted_edges.begin(),sorted_edges.end(),sort_second<igraph_real_t>);

    // Randomize the edges with the same jaccard index. Introduces some variability in
    // the final partitioning.
    size_t i=0;
    while (i<nEdges)
    {
        size_t k=0;
        while ((fabs(sorted_edges.at(i).second - sorted_edges.at(i+k).second)<1E-6) )
        {
            ++k;
            if ((i+k)>=nEdges )
            {
                k-=1;
                break;
            }
        }
        if (k>1)
            std::random_shuffle(sorted_edges.begin()+i,sorted_edges.begin()+(k-1)+i);
        if (k==0)
            break; // no further improvement has been done
        i+=k;
    }
}

/**
 * @brief CommunityStructure::get_sorted_edges_indices
 * @return
 */
vector<int> CommunityStructure::get_sorted_edges_indices()
{
    size_t ms = sorted_edges.size();
    vector<int> ei(ms);
    for (size_t i=0; i<ms; ++i)
        ei.at(i) = sorted_edges.at(i).first;
    return ei;
}

/**
 * @brief CommunityStructure::get_membership
 * @return
 */
const igraph_vector_t* CommunityStructure::get_membership() const
{
    return &this->membership;
}

size_t CommunityStructure::get_membership(size_t i) const
{
    return static_cast<size_t>(std::floor(this->membership.stor_begin[i]));
}

