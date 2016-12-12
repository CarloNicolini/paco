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


#ifndef _IGRAPH_ADDITIONAL_UTILS_
#define _IGRAPH_ADDITIONAL_UTILS_

#include <vector>
#include <igraph.h>
#include <igraph_error.h>
#include <stdexcept>
#include <sstream>

#define IGRAPH_TRY(call){\
    int __result = call;\
    std::stringstream ss; ss << __result; \
    if (__result != 0)\
{\
    throw std::runtime_error("IGRAPH Error " + ss.str());\
    } \
    }

#define IGRAPHPP_TRY_NEW(variable, type)   \
    try {                                \
    variable = new type;             \
    } catch (const std::bad_alloc&) { \
    IGRAPH_ERROR("std::bad_alloc thrown in C++ code", IGRAPH_ENOMEM); \
    }

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols);

igraph_vector_t* order_membership(const igraph_vector_t *curmemb);

int igraph_similarity_jaccard_weighted_pairs(const igraph_t *graph, igraph_vector_t *res,
                                             const igraph_vector_t *pairs, const igraph_vector_t *weights, igraph_neimode_t mode, igraph_bool_t loops);

int igraph_similarity_jaccard_weighted_es(const igraph_t *graph, igraph_vector_t *res,
                                          const igraph_es_t es, const igraph_vector_t *weights, igraph_neimode_t mode, igraph_bool_t loops);

int igraph_i_neisets_intersect(const igraph_t *graph, const igraph_vector_t *v1, const igraph_vector_t *v2, const igraph_vector_t *weights, double *weight_union, double *weight_intersection);

int igraph_read_graph_weighted_edgelist(igraph_t *graph, FILE *instream, igraph_integer_t n, igraph_bool_t directed, igraph_vector_t *edge_weights);

void free_complist(igraph_vector_ptr_t *complist);


#endif
