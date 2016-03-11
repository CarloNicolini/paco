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

#include "igraph_utils.h"

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}

igraph_vector_t* order_membership(const igraph_vector_t *curmemb)
{
        int nVertices = igraph_vector_size(curmemb);
        // groupmap
        igraph_vector_t group_map; igraph_vector_init(&group_map,nVertices);
        igraph_vector_fill(&group_map,-1);
        int current_number=0;
        for (int i=0; i<nVertices; ++i)
        {
            if (group_map.stor_begin[int(curmemb->stor_begin[i])] < 0)
            {
                group_map.stor_begin[int(curmemb->stor_begin[i])] = current_number;
                ++current_number;
            }
        }


        igraph_vector_t *newmemb = new igraph_vector_t;
        igraph_vector_init(newmemb,nVertices);
        igraph_vector_fill(newmemb,-1);

        for (int i=0; i<nVertices; ++i)
        {
            int val = group_map.stor_begin[int(curmemb->stor_begin[i])];
            if (val<0)
                val = *(group_map.stor_end);
            newmemb->stor_begin[i] = val;
        }


        // Clean memory
        igraph_vector_destroy(&group_map);

        return newmemb;
}
