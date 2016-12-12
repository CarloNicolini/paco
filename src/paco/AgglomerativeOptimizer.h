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


#ifndef _AGGLOMERATIVE_OPTIMIZER_H_
#define _AGGLOMERATIVE_OPTIMIZER_H_

#include "QualityOptimizer.h"

class AgglomerativeOptimizer : public QualityOptimizer
{
public:
    AgglomerativeOptimizer() {}
    AgglomerativeOptimizer(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb,const igraph_vector_t *weights=NULL);
    virtual ~AgglomerativeOptimizer();
    double optimize(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, const igraph_vector_t *weights=NULL);
    void set_edges_order(const vector<int> &value);

protected:
    double diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm, const igraph_vector_t *weights);
    vector<int> edges_order;
};

#endif // AgglomerativeOptimizer_H
