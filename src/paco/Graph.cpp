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

#include "Graph.h"

/**
 * @brief GraphC::GraphC
 */
GraphC::GraphC()
{
    IGRAPH_TRY(igraph_empty(&this->ig,0,IGRAPH_UNDIRECTED));
    _is_directed = false;
    _must_delete=true;
}

/**
 * @brief GraphC::GraphC
 * @param g
 */
GraphC::GraphC(igraph_t *g)
{
    IGRAPH_TRY(igraph_empty(&this->ig,0,IGRAPH_UNDIRECTED));
    IGRAPH_TRY(igraph_copy(g,&this->ig));
    _is_directed = false;
    _is_weighted = false;
    _must_delete = true;
}

/**
 * @brief GraphC::GraphC
 * @param nvertices
 */
GraphC::GraphC(size_t nvertices)
{
    IGRAPH_TRY(igraph_empty(&this->ig,nvertices,IGRAPH_UNDIRECTED));
    _is_directed = false;
    _is_weighted = false;
    _must_delete = true;
}

/**
 * @brief GraphC::GraphC
 * @param W
 */
GraphC::GraphC(const Eigen::MatrixXd &W)
{
    this->init(W);
}

/**
 * @brief GraphC::GraphC
 * @param W
 * @param n
 * @param m
 */
GraphC::GraphC(double *W, int n, int m)
{
    Eigen::MatrixXd MW = Eigen::Map<Eigen::MatrixXd>(W,n,m);
    this->init(MW);
}

/**
 * @brief GraphC::init
 * @param W
 */
void GraphC::init(const Eigen::MatrixXd &W)
{
    _must_delete = true;
    igraph_matrix_t w_adj;

    if (W.rows()!=W.cols())
    {
        throw std::logic_error("Non square adjacency matrix");
    }

    if (W.diagonal().sum()>double(0.0))
    {
        //_must_delete = false;
        throw std::logic_error("Adjacency matrix has self-loops, only simple graphs allowed");
    }

    igraph_matrix_view(&w_adj,const_cast<igraph_real_t*>(W.data()),W.cols(),W.rows());
    igraph_weighted_adjacency(&this->ig,&w_adj,IGRAPH_ADJ_UNDIRECTED,NULL,true);

    // Populate the edge weights vector
    edge_weights_stl.clear();
    int n = W.rows();
    int m = W.cols();
    for (int i=0; i<n; ++i)
    {
        for (int j=i+1; j<m; j++)
        {
            double w = W.coeffRef(i*n+j);
            if (w>0)
                edge_weights_stl.push_back(w);
            if (w<0)
            {
                //_must_delete = false;
                throw std::logic_error("Negative edge weight found. Only positive weights supported.");
            }
        }
    }

    if (edge_weights_stl.size() != static_cast<unsigned int>(igraph_ecount(&ig)))
    {
        //_must_delete = false;
        throw std::logic_error("Non consistent length of edge weigths vector, or diagonal entries in adjacency matrix.");
    }

    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size());
    _is_directed = false;

    size_t num_different_edge_weight_values = set<double>(W.data(),W.data()+W.rows()*W.cols()).size();
    _is_weighted = num_different_edge_weight_values!=2;
}

/**
 * @brief GraphC::get_edge_weights
 * @return
 */
const igraph_vector_t* GraphC::get_edge_weights() const
{
    if (_is_weighted)
        return &this->edge_weights;
    else
        throw std::logic_error("Graph is not weighted. Invalid call to get_edge_weights()");
}

/**
 * @brief GraphC::compute_vertex_degrees
 * @param loops
 */
void GraphC::compute_vertex_degrees(bool loops)
{
    igraph_degree(&ig,&vertices_degrees,igraph_vss_all(),IGRAPH_ALL,loops);
}

/**
 * @brief GraphC::compute_vertex_strenghts
 * @param loops
 */
void GraphC::compute_vertex_strenghts(bool loops)
{
    if (is_weighted())
    {
        igraph_strength(&ig,&vertices_strenghts,igraph_vss_all(),IGRAPH_ALL,loops,&edge_weights);
    }
    else
    {
        compute_vertex_degrees(loops);
        cerr << "Vertex strenght = vertex degree on unweighted graph" << endl;
    }
}

/**
 * @brief GraphC::get_degrees
 * @return
 */
const igraph_vector_t* GraphC::get_degrees() const
{
    return &vertices_degrees;
}

/**
 * @brief GraphC::get_strenghts
 * @return
 */
const igraph_vector_t* GraphC::get_strenghts() const
{
    if ( is_weighted())
        return &vertices_strenghts;
    else
        return &vertices_degrees;
}

/**
 * @brief GraphC::~GraphC
 */
GraphC::~GraphC()
{
    if (_must_delete)
        igraph_destroy(&this->ig);
}

/**
 * @brief GraphC::get_igraph
 * @return
 */
const igraph_t* GraphC::get_igraph() const
{
    return &this->ig;
}
#include <sstream>
/**
 * @brief GraphC::read_adj_matrix
 * @param filename
 * @return
 */
bool GraphC::read_adj_matrix(const std::string &filename)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    Eigen::MatrixXd adj_mat;

    ifstream file;
    file.open(filename.c_str(),std::ios::in);

    // Error file is not open
    if(!file.is_open())
    {
        return false;
    }

    vector<vector<double> > data;
    string line;

    while(!std::getline(file, line, '\n').eof())
    {
        std::istringstream reader(line);
        vector<double> lineData;
        while(!reader.eof())
        {
            double val;
            reader >> val;
            if(reader.fail())
                break;
            lineData.push_back(val);
        }
        data.push_back(lineData);
    }
    // Deep copy of data array to the adjacency matrix
    adj_mat.setZero(data.size(),data.size());
    for (unsigned int i=0; i<data.size(); ++i)
        for (unsigned int j=0; j<data.size(); ++j)
            adj_mat.coeffRef(i,j)=data[i][j];

    this->init(adj_mat);
    return true;
}

/**
 * @brief GraphC::read_pajek
 * @param filename
 * @return
 */
bool GraphC::read_pajek(const std::string &filename)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    FILE *f = fopen(filename.c_str(),"r");
    IGRAPH_TRY(igraph_read_graph_pajek(&this->ig,f));
    fclose(f);
    this->_is_weighted = false;
    return true;
}

/**
 * @brief GraphC::read_edge_list
 * @param filename
 * @param nvertices
 * @return
 */
bool GraphC::read_edge_list(const std::string &filename, int nvertices)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    FILE *f = fopen(filename.c_str(),"r");
    IGRAPH_TRY(igraph_read_graph_edgelist(&this->ig,f,nvertices,0));
    fclose(f);
    this->_is_weighted = false;
    return true;
}

/**
 * @brief GraphC::read_weights_from_file
 * @param filename
 * @return
 */
bool GraphC::read_weights_from_file(const string &filename)
{
    ifstream ifs;
    ifs.open(filename.c_str());

    if (!ifs.good())
        throw std::ios_base::failure("Error, file " + filename+ " doesn't exist");

    this->edge_weights_stl.clear();

    string line;
    while ( getline(ifs,line) )
    {
        if (line.size() == 0)
            continue;
        std::stringstream str(line);
        igraph_real_t ew;
        str >> ew;
        if (ew<0)
            throw std::logic_error("Negative edge weight found.");
        edge_weights_stl.push_back(ew);
    }
    ifs.close();

    if (edge_weights_stl.size() != number_of_edges())
        throw std::logic_error("Edge weights vector not consistent with number of graph edges");

    size_t num_different_edge_weight_values = set<double>(edge_weights_stl.begin(),edge_weights_stl.end()).size();
    _is_weighted = num_different_edge_weight_values!=2; // set if the graph is weighted or just has two different weights
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size()); // copy to weights vector
    return true;
}

/**
 * @brief GraphC::set_edge_weights
 * @param w
 */
void GraphC::set_edge_weights(const vector<igraph_real_t> &w)
{
    this->edge_weights_stl = w;
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size());
    size_t num_different_edge_weight_values = set<double>(edge_weights_stl.begin(),edge_weights_stl.end()).size();
    _is_weighted = num_different_edge_weight_values!=2;
}


/**
 * @brief GraphC::read_gml
 * @param filename
 * @return
 */
bool GraphC::read_gml(const string &filename)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    FILE *f = fopen(filename.c_str(),"r");
    IGRAPH_TRY(igraph_read_graph_gml(&this->ig,f));
    fclose(f);
    this->_is_weighted = false;
    return true;
}

/**
 * @brief GraphC::number_of_edges
 * @return
 */
size_t GraphC::number_of_edges() const
{
    return igraph_ecount(&this->ig);
}

/**
 * @brief GraphC::number_of_nodes
 * @return
 */
size_t GraphC::number_of_nodes() const
{
    return igraph_vcount(&this->ig);
}

/**
 * @brief GraphC::add_edge
 * @param source
 * @param target
 * @return
 */
bool GraphC::add_edge(size_t source, size_t target)
{
    if (!is_edge(source,target)) // allow only simple undirected graphs
        igraph_add_edge(&this->ig,source,target);
    return true;
}

/**
 * @brief GraphC::remove_edge
 * @param source
 * @param target
 * @return
 */
bool GraphC::remove_edge(size_t source, size_t target)
{
    if (!is_edge(source,target)) // allow only simple undirected graphs
        return false;
    else
    {
        igraph_es_t es;
        IGRAPH_TRY(igraph_es_pairs_small(&es, IGRAPH_UNDIRECTED, source, target));
        IGRAPH_TRY(igraph_delete_edges(&this->ig,es));
        return true;
    }
}

/**
 * @brief GraphC::remove_edges
 * @param es
 * @return
 */
bool GraphC::remove_edges(igraph_es_t &es)
{
    IGRAPH_TRY(igraph_delete_edges(&this->ig,es));
    return true;
}

/**
 * @brief GraphC::is_edge
 * @param source
 * @param target
 * @return
 */
bool GraphC::is_edge(size_t source, size_t target) const
{
    igraph_bool_t result;
    IGRAPH_TRY(igraph_are_connected(&this->ig, source, target, &result));
    return result;
}

/**
 * @brief GraphC::is_directed
 * @return
 */
bool GraphC::is_directed() const
{
    return igraph_is_directed(&this->ig);
}

/**
 * @brief GraphC::get_edge
 * @param edge_id
 * @return
 */
pair<igraph_integer_t,igraph_integer_t> GraphC::get_edge(size_t edge_id) const
{
    igraph_integer_t from;
    igraph_integer_t to;
    igraph_edge(&this->ig, edge_id,&from,&to);
    return pair<igraph_integer_t,igraph_integer_t>(from,to);
}

/**
 * @brief GraphC::add_vertices
 * @param nvertices
 * @return
 */
bool GraphC::add_vertices(size_t nvertices)
{
    IGRAPH_TRY(igraph_add_vertices(&this->ig,nvertices,0));
    return true;
}

/**
 * @brief GraphC::get_neighbors
 * @param vertex
 * @return
 */
igraph_vector_t GraphC::get_neighbors(size_t vertex)
{
    igraph_vector_t neis;
    IGRAPH_TRY(igraph_neighbors(&this->ig,&neis,vertex,IGRAPH_ALL));
    return neis;
}

/**
 * @brief GraphC::density
 * @param withSelfLoops
 * @return
 */
double GraphC::density(bool withSelfLoops) const
{
    igraph_real_t p;
    IGRAPH_TRY(igraph_density(&this->ig,&p,withSelfLoops));
    return p;
}

/**
 * @brief GraphC::clustering_coefficient
 * @return
 */
double GraphC::clustering_coefficient() const
{
    igraph_real_t clu_coeff;
    IGRAPH_TRY(igraph_transitivity_undirected(&this->ig,&clu_coeff,IGRAPH_TRANSITIVITY_ZERO));
    return clu_coeff;
}


/**
 * @brief GraphC::number_connected_components
 * @return
 */
size_t GraphC::number_connected_components()
{
    long min_num_vertices=1;
    long max_num_components=-1;
    // Need to init components vector pointer, otherwise segfault
    igraph_vector_ptr_t comps;
    igraph_vector_ptr_init(&comps,0);
    IGRAPH_TRY(igraph_decompose(&this->ig, &comps, IGRAPH_WEAK, max_num_components, min_num_vertices));
    long ncomps = igraph_vector_ptr_size(const_cast<igraph_vector_ptr_t*>(&comps));
    igraph_vector_ptr_destroy(&comps);
    return ncomps;
}

/**
 * @brief GraphC::print
 */
void GraphC::print() const
{
    igraph_adjlist_t adjlist;
    igraph_neimode_t mode = IGRAPH_TOTAL;
    igraph_adjlist_init(&this->ig,&adjlist,mode);
    igraph_adjlist_print(&adjlist);
    igraph_adjlist_destroy(&adjlist);
    //igraph_write_graph_dot(&this->ig,stderr);
}

/**
 * @brief GraphC::info
 */
void GraphC::info()
{
    printf(ANSI_COLOR_YELLOW);
    FILE_LOG(logINFO) << "Is directed? " << is_directed();
    FILE_LOG(logINFO) << "Is weighted? " << is_weighted();
    FILE_LOG(logINFO) << "Num vertices=" << number_of_nodes() << " Num edges=" << number_of_edges();
    FILE_LOG(logINFO) << "Density (no loops)=" << density(false);
    FILE_LOG(logINFO) << "Density (with loops)=" << density(true);
    FILE_LOG(logINFO) << "Num connected components=" << number_connected_components();
    FILE_LOG(logINFO) << "Clustering coefficient=" << clustering_coefficient();
    printf(ANSI_COLOR_RESET);
}
