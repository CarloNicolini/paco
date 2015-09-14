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

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <stdexcept>
#include <exception>

#include <Eigen/Core>

#include <igraph.h>
#include "Common.h"
#include "igraph_utils.h"
#include "FileLogger.h"


class GraphC
{
public:
    GraphC();
    GraphC(size_t nvertices);
    GraphC(const Eigen::MatrixXd &W);
    GraphC(double *W, int n, int m);
    GraphC(igraph_t *g);
    ~GraphC();

    void init(const Eigen::MatrixXd &W);

    const igraph_t *get_igraph() const; // allows only const copies of the pointer.
    bool read_pajek(const std::string &filename);
    bool read_edge_list(const std::string &filename, int nvertices);
    bool read_gml(const std::string &filename);
    bool read_weights_from_file(const string &filename);
    void set_edge_weights(const std::vector<igraph_real_t> &w);

    void compute_vertex_strenghts(bool loops=false);
    void compute_vertex_degrees(bool loops=false);
    const igraph_vector_t *get_strenghts() const;
    const igraph_vector_t *get_degrees() const;

    size_t number_of_nodes() const;
    size_t number_of_edges() const;
    double density(bool withSelfLoops=false) const;
    double clustering_coefficient() const;
    igraph_vector_t get_neighbors(size_t vertex);

    bool add_vertices(size_t nvertices);
    bool add_vertices(const std::vector<size_t> &vertices);
    bool add_vertices(igraph_vs_t &vs);
    bool add_edge(size_t source, size_t target);
    bool remove_edge(size_t source, size_t target);
    bool remove_edges(igraph_es_t &es);

    std::pair<igraph_integer_t,igraph_integer_t> get_edge(size_t edgeid) const;
    const igraph_vector_t* get_edge_weights() const;

    bool is_edge(size_t source, size_t target) const;
    bool is_directed() const;
    bool is_weighted() const
    {
        return _is_weighted;
    }
    size_t number_connected_components();

    void print() const;
    void info();

protected:
    igraph_t ig;

    // Weights storage
    vector<igraph_real_t> edge_weights_stl;
    igraph_vector_t edge_weights;

    vector<igraph_real_t> vertices_strenghts_stl;
    igraph_vector_t vertices_strenghts;

    vector<igraph_real_t> vertices_degrees_stl;
    igraph_vector_t vertices_degrees;

    bool _is_weighted;
    bool _is_directed;
};

#endif
