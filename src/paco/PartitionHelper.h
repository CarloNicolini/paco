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

#ifndef PARTITIONHELPER_H
#define PARTITIONHELPER_H

#include <algorithm> // for std::count
#include <igraph.h>

#include "igraph_utils.h"
#include "Common.h"

typedef unordered_map<size_t,size_t> IntMap;
typedef unordered_map<size_t,double> DoubleMap;
typedef unordered_map<size_t, set<size_t> > CommMap;
typedef unordered_map<size_t, set<size_t> >::iterator CommMapIter;
typedef unordered_map<size_t, set<size_t> >::const_iterator CommMapCIter;

class PartitionHelper
{
public:
    PartitionHelper();
    ~PartitionHelper();

    void init(const igraph_t*g, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    bool move_vertex(const igraph_t *g, const igraph_vector_t *memb, int vert, size_t dest_comm, const igraph_vector_t *weights=NULL);
    bool merge_communities(const igraph_t *g, const igraph_vector_t *memb, size_t source_comm, size_t dest_comm, const igraph_vector_t *weights=NULL);
    inline size_t get_membership(const igraph_vector_t *memb, int vert) const;
    void reindex(const igraph_vector_t *memb);
    void print() const;

    const double& get_graph_total_pairs() const
    {
        return graph_total_pairs;
    }

    const double& get_graph_total_weight() const
    {
        return graph_total_weight;
    }

    const double& get_total_incomm_pairs() const
    {
        return total_incomm_pairs;
    }

    const double& get_total_incomm_weight() const
    {
        return total_incomm_weight;
    }

    const size_t& get_num_vertices() const
    {
        return num_vertices;
    }

    const size_t& get_num_edges() const
    {
        return num_edges;
    }

    const igraph_vector_t* get_all_degrees() const
    {
        return &all_degrees;
    }

    size_t get_num_comms() const
    {
        return num_comms;
    }

    const IntMap &get_incomm_nvert() const
    {
        return incomm_nvert;
    }

    const DoubleMap &get_incomm_weight() const
    {
        return incomm_weight;
    }

    const IntMap &get_incomm_pairs() const
    {
        return incomm_pairs;
    }

    const DoubleMap &get_incomm_deg() const
    {
        return incomm_deg;
    }

    const CommMap& get_communities() const
    {
        return communities;
    }

protected:
    igraph_vector_t all_degrees;
    igraph_vector_t all_strenght;

    IntMap incomm_nvert;
    IntMap incomm_pairs;
    DoubleMap incomm_weight;
    DoubleMap incomm_deg;

    const igraph_vector_t *curmemb;

    CommMap communities;
    size_t num_comms;
    size_t num_vertices;    // number of vertices
    size_t num_edges; // number of edges
    double total_incomm_weight; // sum of all edge weights inside communities
    double total_incomm_pairs;  // sum of all vertex pairs inside communities

    double graph_total_weight; // number of edges
    double graph_total_pairs; // // n*(n-1)/2, number of total graph vertex pairs

private:
    inline bool check_comm(int dest_comm);
    void fill_communities(const igraph_vector_t *memb);
    double weight_to_from_community(const igraph_t *g, const igraph_vector_t* memb, size_t v, size_t comm, igraph_neimode_t mode, const igraph_vector_t *weights=NULL);
};


#endif // PARTITIONHELPER_H
