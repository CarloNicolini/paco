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

#include "QualityFunction.h"
#include "DegreeCorrectedSurpriseFunction.h"
#include "igraph_utils.h"
#include "KLDivergence.h"

/**
 * @brief DegreeCorrectedSurpriseFunction::DegreeCorrectedSurpriseFunction
 */
DegreeCorrectedSurpriseFunction::DegreeCorrectedSurpriseFunction() {}

/**
 * @brief DegreeCorrectedSurpriseFunction::eval
 * @param g
 * @param memb
 * @param weights
 */
void DegreeCorrectedSurpriseFunction::eval(const igraph_t *g, const igraph_vector_t *memb, const igraph_vector_t *weights) const
{
    PartitionHelper *par = new PartitionHelper;
    par->init(g,memb,weights);
    this->eval(par);
    delete par;
    //printf("--> AS=%f -- mi=%f pi=%f m=%f p=%f\n",quality, mi,pi,m,p);
}

/**
 * @brief DegreeCorrectedSurpriseFunction::eval
 * @param par
 */
void DegreeCorrectedSurpriseFunction::eval(const PartitionHelper *par) const
{
    quality = 0;
    igraph_real_t m = par->get_graph_total_weight();
    igraph_real_t p = par->get_graph_total_pairs();

    for (CommMapCIter iter=par->get_communities().begin(); iter!=par->get_communities().end(); ++iter)
    {
        size_t c = iter->first;
        igraph_real_t mc = par->get_incomm_weight().at(c);
        igraph_real_t Kc = par->get_incomm_deg().at(c);
        igraph_real_t nc = par->get_incomm_nvert().at(c);
        igraph_real_t pc = nc*(nc-1.0)/2.0;
        for (CommMapCIter iter=par->get_communities().begin(); iter!=par->get_communities().end(); ++iter)
        {
            size_t d = iter->first;
            igraph_real_t mcd = 0;//par->weight_to_from_community(par->get_i)
            igraph_real_t Kd = par->get_incomm_deg().at(d);
            quality *= logC(Kc*Kd,mcd)-logC(4*m*m,2*m);
        }
    }
}
