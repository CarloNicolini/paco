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


#include <sstream>
#include "Graph.h"
#include "FileLogger.h"
/**
 * @brief GraphC::GraphC
 */
GraphC::GraphC()
{
    IGRAPH_TRY(igraph_empty(&this->ig,0,IGRAPH_UNDIRECTED));
    _is_directed = false;
    _is_weighted = false;
    _has_selfloops = false;
    _must_delete=true;
}

/**
 * @brief GraphC::GraphC
 * @param rhs
 */
GraphC::GraphC(const GraphC &rhs)
{
    igraph_copy(&this->ig,rhs.get_igraph());

    // Copy other private internals
    this->_is_weighted = rhs._is_weighted;
    this->_is_directed = rhs._is_directed;
    this->_must_delete = rhs._must_delete;
    this->_has_selfloops = rhs._has_selfloops;

    // Copy edge weights
    this->edge_weights_stl = rhs.edge_weights_stl;
    // Copy vertices strenghts
    this->vertices_strenghts_stl = rhs.vertices_strenghts_stl;
    // Copy vertices degrees
    this->vertices_degrees_stl = rhs.vertices_degrees_stl;

    // Set the edge weights vector as a view on underlying STL structure
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size());
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
 * @brief GraphC::GraphC
 * @param edges_list_stl
 * @param weights
 */
GraphC::GraphC(const double *edges_list, const double *edges_weights, int nedges)
{
    this->init(edges_list,edges_weights,nedges);
}

/**
 * @brief GraphC::init
 * @param ewlist
 * @param num_edges
 */
void GraphC::init(const double *ewlist, int _weighted, int num_edges)
{
    vector<double> elist,wlist;
    int iplus = (_weighted? 3 : 2);
    for (int i=0; i<iplus*num_edges; i+=iplus)
    {
        elist.push_back(ewlist[i]);
        elist.push_back(ewlist[i+1]);
    }
    for (int i=0; i<iplus*num_edges; i+=iplus)
    {
        if (_weighted)
            wlist.push_back(ewlist[i+2]);
        else
            wlist.push_back(1);
    }
    this->init(elist.data(),wlist.data(),num_edges);
}

/**
 * @brief GraphC::init
 * @param elist
 * @param weights
 * @param num_edges
 */
void GraphC::init(const double *elist, const double *weights, int num_edges)
{
    igraph_vector_t edges_list;
    igraph_vector_view(&edges_list,elist,2*num_edges);
    igraph_create(&this->ig, &edges_list, 0, 0);

    // Assing weights
    vector<igraph_real_t> weights_stl(num_edges,1.0);
    if (weights!=NULL)
        weights_stl.assign(weights, weights + num_edges);
    
    this->set_edge_weights(weights_stl);
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size());
    _must_delete = true;
    // Check for self-loops
    _is_directed = this->ig.directed;
    for (int i=0; i<num_edges*2; i+=2)
    {
        _has_selfloops = (elist[i]==elist[i+1]);
    }
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
    _is_weighted = num_different_edge_weight_values != 2;
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
        return NULL;
}

/**
 * @brief GraphC::compute_vertex_degrees
 * @param loops
 */
void GraphC::compute_vertex_degrees(bool loops)
{
    if (this->number_of_nodes()==0)
        return;
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
    if (igraph_vector_size(&vertices_degrees)==0)
    {
        throw std::runtime_error("Vertex degrees vector is empty. Call compute_vertex_degrees() first");
    }
    return &vertices_degrees;
}

/**
 * @brief GraphC::get_strenghts
 * @return
 */
const igraph_vector_t* GraphC::get_strenghts() const
{
    if ( is_weighted())
    {
        if (igraph_vector_size(&vertices_strenghts)==0)
        {
            throw std::runtime_error("Vertex strenghts vector is empty. Call compute_vertex_strenghts() first");
        }
        return &vertices_strenghts;
    }
    else
    {
        if (igraph_vector_size(&vertices_degrees)==0)
        {
            throw std::runtime_error("Vertex degrees vector is empty. Call compute_vertex_degrees() first");
        }
        return &vertices_degrees;
    }
}

/**
 * @brief GraphC::~GraphC
 */
GraphC::~GraphC()
{
    if (_must_delete)
    {
        igraph_destroy(&this->ig);
    }
}

/**
 * @brief GraphC::get_igraph
 * @return
 */
const igraph_t* GraphC::get_igraph() const
{
    return &this->ig;
}

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
bool GraphC::read_edge_list(const std::string &filename)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    FILE *f = fopen(filename.c_str(),"r");
    IGRAPH_TRY(igraph_read_graph_edgelist(&this->ig,f,0,0));
    fclose(f);
    this->_is_weighted = false;
    return true;
}

/**
 * @brief GraphC::read_weighted_edge_list
 * @param filename
 * @param nvertices
 * @param edge_weights vector
 * @return
 */
bool GraphC::read_weighted_edge_list(const std::string &filename)
{
    IGRAPH_TRY(igraph_destroy(&this->ig));
    FILE *f = fopen(filename.c_str(),"r");
    IGRAPH_TRY(igraph_read_graph_weighted_edgelist(&this->ig,f,0,0,&this->edge_weights));
    fclose(f);
    this->_is_weighted = true;
    int m = igraph_vector_size(&this->edge_weights);
    this->edge_weights_stl.resize(m);
    for (int i=0; i<m; ++i)
        edge_weights_stl.at(i) = edge_weights.stor_begin[i];
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
    _is_weighted = num_different_edge_weight_values != 2; // set if the graph is weighted or just has two different weights
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size()); // copy to weights vector
    return true;
}

/**
 * @brief GraphC::set_edge_weights
 * @param w
 */
void GraphC::set_edge_weights(const vector<igraph_real_t> &w, bool override_is_weighted)
{
    this->edge_weights_stl = w;
    igraph_vector_view(&edge_weights,edge_weights_stl.data(),edge_weights_stl.size());

    size_t num_different_edge_weight_values = set<double>(edge_weights_stl.begin(),edge_weights_stl.end()).size();
    _is_weighted = num_different_edge_weight_values > 1;
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
 * @brief GraphC::line_graph Create a line graph from the current graph.
 * The line graph L(G) of a G undirected graph is defined as follows. L(G) has one vertex
 * for each edge in G and two vertices in L(G) are connected by an edge if their
 * corresponding edges share an end point.
 * @return a pointer to a GraphC representing the line graph. Must be deleted.
 */
GraphC* GraphC::line_graph()
{
    igraph_t *lg=NULL;
    igraph_linegraph(&(this->ig),lg);
    return new GraphC(lg);
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

    igraph_vector_ptr_t complist;
    igraph_vector_ptr_init(&complist, 0);
    igraph_decompose(&this->ig, &complist, IGRAPH_WEAK, max_num_components, min_num_vertices);
    long ncomps = igraph_vector_ptr_size(&complist);
    free_complist(&complist);
    igraph_vector_ptr_destroy(&complist);
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
    FILE_LOG(logINFO) << "Graph adjacencly list" ;
    igraph_adjlist_print(&adjlist);
    FILE_LOG(logINFO) << "End adjacencly list" ;
    igraph_adjlist_destroy(&adjlist);
}

/**
 * @brief GraphC::info
 */
void GraphC::info()
{
    printf(ANSI_COLOR_YELLOW);
    FILE_LOG(logINFO) << "Num vertices=" << number_of_nodes() << " Num edges=" << number_of_edges();
    FILE_LOG(logINFO) << "Is directed? " << is_directed();
    FILE_LOG(logINFO) << "Is weighted? " << is_weighted();
    FILE_LOG(logINFO) << "Density (no loops)=" << density(false);
    FILE_LOG(logINFO) << "Density (with loops)=" << density(true);
    FILE_LOG(logINFO) << "Num connected components=" << number_connected_components();
    FILE_LOG(logINFO) << "Clustering coefficient=" << clustering_coefficient();
    printf(ANSI_COLOR_RESET);
}
