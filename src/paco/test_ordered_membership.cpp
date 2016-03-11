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

#include <iostream>

#include "Graph.h"
#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "ModularityFunction.h"
#include "Community.h"
#include "PartitionHelper.h"
#include "igraph_utils.h"

using namespace std;

int main(int argc, char *argv[])
{
/*
    GraphC h;
    h.read_adj_matrix(string(argv[1]));
    const igraph_t *g = h.get_igraph();
    CommunityStructure c(&h);
    c.read_membership_from_file(argv[2]);

    const igraph_vector_t *m = c.get_membership();

    c.order_membership();
*/
    double v[34] = {4, 7, 3, 3, 4, 4, 4, 3, 8, 9, 4, 11, 12, 3, 8, 8, 4, 17, 8, 19, 8, 21, 22, 27, 31, 25, 26, 27, 25, 26, 8, 25, 8, 8};

    const igraph_vector_t *memb = new igraph_vector_t;
    memb = igraph_vector_view(memb,v,34);

    igraph_vector_t *newmemb = order_membership(memb);
    igraph_vector_print(newmemb);
    return 0;
}
