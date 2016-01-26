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
#include "SignificanceFunction.h"
#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "Community.h"
#include "PartitionHelper.h"
#include "RandomOptimizer.h"
#include "AnnealOptimizer.h"
#include "AgglomerativeOptimizer.h"
#include "Timer.h"
#include "ModularityFunction.h"
#include "AnnealOptimizer.h"
using namespace std;

int main(int argc, char *argv[])
{
    GraphC h;
    h.read_gml(string(argv[1]));

    CommunityStructure c(&h);
    c.set_random_seed();
    c.sort_edges();

    AgglomerativeOptimizer opt;

    AsymptoticSurpriseFunction fun;
    //c.sort_edges();
    //opt.set_edges_order(c.get_sorted_edges_indices());

    //for (int i=0; i<1;i++)
    //    cerr << "AS=" << opt.optimize(h.get_igraph(),fun,c.get_membership(),NULL) << endl;

    AnnealOptimizer ann;
    AnnealParameters pars;
    pars.temperature = 1;
    ann.set_parameters(pars);

    cerr << "Annealing AS=" << ann.optimize(h.get_igraph(),fun,c.get_membership(),NULL) << endl;

    c.reindex_membership();
    c.save_membership("test.csv");

    igraph_vector_print(c.get_membership());

    return 0;
}
