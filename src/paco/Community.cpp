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


#include <igraph.h>
#include "Community.h"

#include "QualityFunction.h"
#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "SignificanceFunction.h"

#ifdef EXPERIMENTAL
#include "AsymptoticModularityFunction.h"
#include "WonderFunction.h"
#include "DegreeCorrectedSurpriseFunction.h"
#endif

#include "QualityOptimizer.h"
#include "AnnealOptimizer.h"
#include "AgglomerativeOptimizer.h"
#include "RandomOptimizer.h"

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
    if (vmemb.size() != static_cast<unsigned int>(this->nVertices))
        throw std::runtime_error("Error loading membership file. Vertex number not consistent with current graph.");

    for (igraph_integer_t i=0; i<this->nVertices; ++i)
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
        for (long int i=0; i<igraph_vector_size(m); i++)
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
    if (G->number_of_nodes() ==0)
        throw std::logic_error("Graph with no vertices");
    if (G->number_of_edges()==0)
        throw std::logic_error("Graph with no edges");
    this->nVertices = G->number_of_nodes();
    this->nEdges = G->number_of_edges();
    this->nPairs = nVertices*(nVertices-1)/2;

    // Init the membership vector
    IGRAPH_TRY(igraph_vector_init(&membership,nVertices));
    for (igraph_integer_t i=0; i<nVertices; ++i)
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
    for (igraph_integer_t i=0; i<nVertices; ++i)
        VECTOR(membership)[i] -= minC;
}

/**
 * @brief CommunityStructure::order_membership
 * Order the membership vector such that
 * c_i <= max(c_1,c_2,...,c_(i-i))+1 (see http://tuvalu.santafe.edu/~aaronc/modularity/ for details)
 */
void CommunityStructure::order_membership()
{
    // groupmap
    igraph_vector_t group_map; igraph_vector_init(&group_map,this->nVertices);
    igraph_vector_fill(&group_map,-1);
    int current_number=0;
    for (int i=0; i<this->nVertices; ++i)
    {
        if (group_map.stor_begin[int(membership.stor_begin[i])] < 0)
        {
            group_map.stor_begin[int(membership.stor_begin[i])] = current_number;
            ++current_number;
        }
    }


    igraph_vector_t newmemb;
    igraph_vector_init(&newmemb,this->nVertices);
    igraph_vector_fill(&newmemb,-1);

    for (int i=0; i<this->nVertices; ++i)
    {
        int val = group_map.stor_begin[int(membership.stor_begin[i])];
        if (val<0)
            val = *(group_map.stor_end);
        newmemb.stor_begin[i] = val;
    }

    //igraph_vector_print(&group_map);
    //igraph_vector_print(&newmemb);

    igraph_vector_update(&membership,&newmemb);
    // Clean memory
    igraph_vector_destroy(&group_map);
    igraph_vector_destroy(&newmemb);
}


/**
 * @brief CommunityStructure::optimize
 * @param qual
 * @param optmethod
 * @param nrep
 * @return
 */
double CommunityStructure::optimize(QualityType qual, OptimizerType optmethod, int nrep)
{
    QualityOptimizer *opt;
    QualityFunction *fun;

    const igraph_vector_t *edge_weights = pgraph->get_edge_weights();

    double finalqual = std::numeric_limits<double>::min();

    // Select the quality function
    switch (qual)
    {
    case QualitySurprise:
    {
        if (pgraph->is_weighted())
        {
            throw std::logic_error("Can't optimize discrete surprise on weighted graph. Use AsymptoticSurprise instead.");
        }
        else
        {
            fun = dynamic_cast<SurpriseFunction*>(new SurpriseFunction);
        }
        break;
    }
    case QualitySignificance:
    {
        fun = dynamic_cast<SignificanceFunction*>(new SignificanceFunction);
        break;
    }
    case QualityAsymptoticSurprise:
    {
        fun = dynamic_cast<AsymptoticSurpriseFunction*>(new AsymptoticSurpriseFunction);
        break;
    }
    case QualityInfoMap:
    {
        // For Infomap it selects the partition with the minimum description length (last argument) and it saves it to final qual
        igraph_community_infomap(pgraph->get_igraph(),edge_weights,NULL,nrep,&membership,&finalqual);
        return finalqual;
    }
//    case QualityAsymptoticModularity:
//    {
//        fun = dynamic_cast<AsymptoticModularityFunction*>(new AsymptoticModularityFunction);
//        break;
//    }
//    case QualityWonder:
//    {
//        fun = dynamic_cast<WonderFunction*>(new WonderFunction);
//        break;
//    }
//    case QualityDegreeCorrectedSurprise:
//    {
//        fun = dynamic_cast<DegreeCorrectedSurpriseFunction*>(new DegreeCorrectedSurpriseFunction);
//        break;
//    }
    default:
    {
        throw std::logic_error("Non supported quality function");
    }
    }

    // Select the optimization method
    switch (optmethod)
    {
    case MethodAgglomerative:
    {
        opt = dynamic_cast<AgglomerativeOptimizer*>( new AgglomerativeOptimizer);
        dynamic_cast<AgglomerativeOptimizer*>(opt)->set_edges_order(this->get_sorted_edges_indices());
        break;
    }
    case MethodRandom:
    {
        opt = dynamic_cast<RandomOptimizer*>(new RandomOptimizer);
        break;
    }
    case MethodAnneal:
    {
        opt = dynamic_cast<AnnealOptimizer*>(new AnnealOptimizer);
        break;
    }
    default:
    {
        throw std::logic_error("Non supported optimization method");
    }
    }

    // Now select the partition with the MAXIMUM quality value
    igraph_vector_t best_membership;
    igraph_vector_init(&best_membership,pgraph->number_of_nodes());
    //try
    //{
        for (int i=0; i<nrep; ++i)
        {
            double qual = opt->optimize(pgraph->get_igraph(),*fun,&membership,edge_weights);
            if (qual>finalqual)
            {
                finalqual = qual;
                igraph_vector_update(&best_membership,&membership);
            }
        }
        // then copy back the content of best_membership to membership
        igraph_vector_update(&membership,&best_membership);
        igraph_vector_destroy(&best_membership);
    //}
    //catch ( std::exception &e )
    //{
    //    throw e;
    //}

    delete opt;
    delete fun;
    return finalqual;
}

/**
 * @brief CommunityStructure::compute_edges_similarities
 * See the implementation of "Density-based shrinkage for revealing hierarchical and overlapping
community structure in networks" Physica A 390 (2011) 2160-2171
 */
void CommunityStructure::compute_edges_similarities()
{
    if (!this->pgraph->is_weighted())
    {
        // Query the similarities for all edges
        igraph_similarity_jaccard_es(this->pgraph->get_igraph(),
                                     &edges_sim,
                                     igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                     IGRAPH_ALL,
                                     false);
    }
    else
    {
        igraph_similarity_jaccard_weighted_es(this->pgraph->get_igraph(),&edges_sim,
                                              igraph_ess_all(IGRAPH_EDGEORDER_ID),
                                              this->pgraph->get_edge_weights(),
                                              IGRAPH_ALL,
                                              false);
        /*
        cerr << "=== EDGE WEIGHTS ===" << endl;
        for (igraph_integer_t i=0; i<this->nEdges; ++i)
        {
            cout << "(" << this->pgraph->get_edge(i).first << " " << this->pgraph->get_edge(i).second << ")" ;
            cout << "\t" << this->pgraph->get_edge_weights()->stor_begin[i] << endl;
        }
        */
    }
}

/**
 * @brief CommunityStructure::sort_edges
 */
void CommunityStructure::sort_edges()
{
    this->compute_edges_similarities(); // initialize the edges_sim igraph_vector_t
    sorted_edges.resize(nEdges);

    for (igraph_integer_t i=0; i<this->nEdges; ++i)
    {
        sorted_edges.at(i).first = i;
        sorted_edges.at(i).second = edges_sim.stor_begin[i];
        //if (this->pgraph->is_weighted()) //XXX Jaccard index is then weighted by edge weight
       //     sorted_edges.at(i).second  *= this->pgraph->get_edge_weights()->stor_begin[i];
    }

    // Edges are sorted by jaccard index. Edges endpoints can be found by IGRAPH_TO and IGRAPH_FROM macros
    std::sort(sorted_edges.begin(),sorted_edges.end(),sort_second<igraph_real_t>);

    /*
    cerr << "=== JACCARD WEIGHTED ===" << endl;
    for (igraph_integer_t i=0; i<nEdges; ++i)
        cerr << "(" << pgraph->get_edge(sorted_edges.at(i).first).first << "," << pgraph->get_edge(sorted_edges.at(i).first).second << ") jacc_sim=" << edges_sim.stor_begin[i] << endl;
    */
    // Randomize the edges with the same jaccard index. Introduces some variability in
    // the final partitioning.
    igraph_integer_t i=0;
    while (i<nEdges)
    {
        igraph_integer_t k=0;
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

/**
 * @brief CommunityStructure::get_membership_vector
 * @return
 */
vector<int> CommunityStructure::get_membership_vector() const
{
    vector<int> memb(this->nVertices);
    for (int i=0; i<nVertices; ++i)
        memb.at(i) = this->membership.stor_begin[i];
    return memb;
}

/**
 * @brief CommunityStructure::get_membership
 * @param i
 * @return
 */
size_t CommunityStructure::get_membership(size_t i) const
{
    return static_cast<size_t>(std::floor(this->membership.stor_begin[i]));
}

