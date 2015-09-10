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


#include "PartitionHelper.h"
#include "igraph_utils.h"

/**
 * @brief PartitionHelper::PartitionHelper
 */
PartitionHelper::PartitionHelper()
{
    igraph_vector_init(&all_degrees,0);
    igraph_vector_init(&all_strenght,0);

    graph_total_weight = 0;
    total_incomm_weight = 0;
    total_incomm_pairs = 0;
    num_vertices = 0;
    graph_total_pairs = 0;
    num_comms = 0;

    curmemb = NULL;
}

/**
 * @brief PartitionHelper::~PartitionHelper
 */
PartitionHelper::~PartitionHelper()
{
    igraph_vector_destroy(&all_degrees);
    igraph_vector_destroy(&all_strenght);
}

/**
 * @brief PartitionHelper::init
 * @param graph
 * @param reind_memb
 * @param weights
 */
void PartitionHelper::init(const igraph_t*graph, const igraph_vector_t *memb, const igraph_vector_t *weights)
{
    this->communities.clear();
    this->incomm_deg.clear();
    this->incomm_nvert.clear();
    this->incomm_pairs.clear();
    this->incomm_weight.clear();
    this->graph_total_pairs = 0;
    this->graph_total_weight = 0;
    this->num_comms = 0;
    this->total_incomm_pairs = 0;
    this->total_incomm_weight = 0;

    this->curmemb = memb;

    this->num_edges = igraph_ecount(graph);
    igraph_real_t m=num_edges;

    this->num_vertices = igraph_vcount(graph);
    this->graph_total_pairs = num_vertices*(num_vertices-1)/2;
    this->graph_total_weight = m;

    if ( igraph_vector_size(memb) < num_vertices )
    {
        throw std::logic_error("Cannot calculate modularity, inconsistent membership vector length");
    }

    igraph_vector_resize(&all_degrees,num_vertices);
    igraph_degree(graph,&all_degrees,igraph_vss_all(),IGRAPH_TOTAL,false);
    // Compute the vector strenghts
    igraph_vector_resize(&all_strenght,num_vertices);
    igraph_strength(graph,&all_strenght,igraph_vss_all(),IGRAPH_TOTAL,false,weights);

    this->fill_communities(memb);

    // Initialize STL map containers with zeros except for sincomm_nvertices,sincomm_pairs with number of vertices/pairs
    for (CommMapCIter it = communities.cbegin(); it!=communities.cend(); ++it)
    {
        size_t c = it->first;
        incomm_weight[c]=0;
        incomm_deg[c]=0;
        incomm_nvert[c]=it->second.size();
        incomm_pairs[c]=num_pairs(it->second.size());
    }

    // Fill sincomm_weight and sincomm_degree
    igraph_integer_t from=0;
    igraph_integer_t to=0;
    size_t c1=0;
    size_t c2=0;

    if (weights)
    {
        if (igraph_vector_size(weights) != m)
            throw std::logic_error("Weights vector != number of edges");
        graph_total_weight = igraph_vector_sum(weights);
        for (size_t ei=0; ei<m; ei++)
        {
            igraph_real_t w= weights->stor_begin[ei];
            if (w < 0)
                throw std::logic_error("Negative weight in weight vector");
            igraph_edge(graph, (igraph_integer_t) ei, &from, &to);
            c1=(size_t) memb->stor_begin[from];
            c2=(size_t) memb->stor_begin[to];
            if (c1==c2)
            {
                incomm_weight[c1]+=w;
                total_incomm_weight += w;
            }
            incomm_deg[c1] += w;
            incomm_deg[c2] += w;
        }
    }
    else
    {
        for (size_t ei=0; ei<m; ei++)
        {
            igraph_edge(graph, (igraph_integer_t) ei, &from, &to);
            c1=(size_t) memb->stor_begin[from];
            c2=(size_t) memb->stor_begin[to];
            if (c1==c2)
            {
                incomm_weight[c1]+=1.0;
                total_incomm_weight += 1.0;
            }
            incomm_deg[c1] += 1.0;
            incomm_deg[c2] += 1.0;
        }
    }

    // Sum total intracluster pairs
    total_incomm_pairs = mapvalue_sum<size_t>(incomm_pairs);
    //total_incomm_weight = mapvalue_sum<double>(incomm_weight); //already computed
}

/**
 * @brief PartitionHelper::count_communities
 * @param memb
 * @return
 */
void PartitionHelper::fill_communities(const igraph_vector_t *memb)
{
    size_t mlen = igraph_vector_size(memb);
    for (size_t c=0; c < mlen; c++)
        communities[memb->stor_begin[c]].insert(c); // insert node c into community memb[c]

    this->num_comms = communities.size();
}

/**
 * @brief PartitionHelper::reindex
 * @param memb
 */
void PartitionHelper::reindex(const igraph_vector_t *memb)
{
    IGRAPH_TRY(igraph_reindex_membership(const_cast<igraph_vector_t *>(memb),NULL));
    int minC = igraph_vector_min(memb)-1; // so to start from 1 to |C| included
    for (size_t i=0; i<igraph_vector_size(memb); ++i)
        memb->stor_begin[i] -= minC;
}

/**
 * @brief PartitionHelper::get_membership
 * @param memb
 * @param vert
 * @return
 */
inline size_t PartitionHelper::get_membership(const igraph_vector_t *memb, int vert) const
{
    if (vert >= num_vertices)
        throw std::range_error("Error indexing vertex");
    return memb->stor_begin[vert];
}

/**
 * @brief PartitionHelper::check_comm
 * @param dest_comm
 * @return
 */
inline bool PartitionHelper::check_comm(int dest_comm)
{
    if (communities.count(dest_comm) == 0)
        return false;
    else
        return true;
}

/**
 * @brief PartitionHelper::move_vertex Move a vertex source to a dest_comm community
 * @param g
 * @param memb
 * @param source
 * @param dest_comm
 * @param weights
 * @return
 */
bool PartitionHelper::move_vertex(const igraph_t *g, const igraph_vector_t * memb, int source, size_t dest_comm, const igraph_vector_t *weights)
{
    size_t source_comm = get_membership(memb,source);
    if (source_comm==dest_comm)
        return false; // do nothing because same community

    if (!check_comm(dest_comm))
    {
        // if dest_comm does not exist, then create an empty dest_comm
        communities[dest_comm] = std::set<size_t>();
        incomm_deg[dest_comm] = 0;
        incomm_nvert[dest_comm] = 0;
        incomm_weight[dest_comm] = 0;
        incomm_pairs[dest_comm] = 0;
    }

    // Update community vector, remove vertex "source" from its original community
    if (communities.at(source_comm).find(source) != communities.at(source_comm).end())
        communities.at(source_comm).erase(source);
    else // if it's not in its original community then throw error
        throw std::logic_error("Vertex not found in source community");
    // and then add it to the destination community if it is not
    if (communities.at(dest_comm).find(source) == communities.at(dest_comm).end() )
        communities.at(dest_comm).insert(source);
    else // if it's already in the destination community, then throw error
        throw std::logic_error("Vertex already in destination community");

    // Compute the number of neighbors that vertex source has in its original community
    double w_in = weight_to_from_community(g,memb,source,source_comm,IGRAPH_ALL,weights);
    // Compute the number of neighbors that vertex source has in destination community
    double w_to = weight_to_from_community(g,memb,source,dest_comm,IGRAPH_ALL,weights);

    incomm_weight.at(source_comm) -= w_in;
    incomm_weight.at(dest_comm) += w_to;

    // Update community num_vertices after movement of vertex source to dest_comm
    size_t old_n_source = incomm_nvert.at(source_comm);
    size_t old_n_dest = incomm_nvert.at(dest_comm);
    int old_p_source = num_pairs(old_n_source);
    int old_p_dest = num_pairs(old_n_dest); // important to use int otherwise values are casted and -1 in pairs doesn't work

    incomm_nvert.at(source_comm) -= 1;//communities.at(source_comm).size(); // this is faster that set.size()
    incomm_nvert.at(dest_comm) +=1;// communities.at(dest_comm).size();

    // Update number of pairs in communities
    incomm_pairs.at(source_comm) = num_pairs(incomm_nvert.at(source_comm));
    incomm_pairs.at(dest_comm) = num_pairs(incomm_nvert.at(dest_comm));

    // Update total intracluster weight
    this->total_incomm_weight += w_to - w_in;

    // Update total intracluster pairs, by computing delta_num_pairs for move_vertex
    int deltap = ( incomm_pairs.at(dest_comm)-old_p_source) + (incomm_pairs.at(source_comm)-old_p_dest);
    total_incomm_pairs += deltap;

    // Finally do the movement!
    memb->stor_begin[source] = dest_comm;
    // Assign current membership pointer
    this->curmemb = memb;
    return true;
}

bool PartitionHelper::merge_communities(const igraph_t *g, const igraph_vector_t *memb, size_t source_comm, size_t dest_comm, const igraph_vector_t *weights)
{
    if (source_comm==dest_comm)
        return false; // do nothing because same community

    if (!check_comm(dest_comm))
        throw std::runtime_error("Non existing destination community");

    set<size_t>::const_iterator it = communities.at(source_comm).cbegin();
    while (it!=communities.at(source_comm).end())
    {
        size_t v = *(it++);
        move_vertex(g,memb,v,dest_comm,weights);
    }
}

/**
 * @brief PartitionHelper::weight_to_from_community
 * @param g
 * @param memb
 * @param v
 * @param comm
 * @param mode
 * @param weights
 * @return
 */
double PartitionHelper::weight_to_from_community(const igraph_t *g, const igraph_vector_t *memb, size_t v, size_t comm, igraph_neimode_t mode, const igraph_vector_t *weights)
{
    if (weights)
    {
        if (igraph_vector_size(weights) != igraph_ecount(g))
            throw std::logic_error("Weights vector != number of edges");
    }

    double total_w = 0.0;
    size_t degree = all_degrees.stor_begin[v];
    igraph_vector_t incident_edges;
    igraph_vector_t neighbours;
    igraph_vector_init(&incident_edges, degree);
    igraph_vector_init(&neighbours, degree);
    igraph_incident(g, &incident_edges, v, mode);
    igraph_neighbors(g, &neighbours, v, mode);

    for (size_t i = 0; i < degree; i++)
    {
        size_t u = VECTOR(neighbours)[i];
        // If it is an edge to the requested community
        if ( get_membership(memb,u) == comm)
        {
            size_t e = VECTOR(incident_edges)[i];
            if (weights)
                total_w += weights->stor_begin[e];
            else
                total_w += 1;
        }
    }
    igraph_vector_destroy(&incident_edges);
    igraph_vector_destroy(&neighbours);

    return total_w;
}

/**
 * @brief PartitionHelper::print
 */
void PartitionHelper::print() const
{
    if (curmemb==NULL)
        throw std::runtime_error("Non initialized membership vector. Call PartitionHelper::init first");

    printf(ANSI_COLOR_BLUE);
    printf("memb:");
    igraph_vector_print(curmemb);

    printf(ANSI_COLOR_YELLOW);
    for (CommMapCIter it = communities.cbegin(); it!=communities.cend(); ++it)
    {
        size_t c = it->first;
        if (it->second.empty())
            continue;

        printf("%zu\t%g\t%zu\t%zu\t{",c,incomm_weight.at(c),incomm_nvert.at(c),incomm_pairs.at(c));
        for ( std::set<size_t>::const_iterator it2 = communities.at(c).begin(); it2!=communities.at(c).end(); ++it2)
        {
            printf("%zu,",*it2);
        }
        printf("}\n");
    }
    printf(ANSI_COLOR_RESET);
}

