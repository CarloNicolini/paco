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

using namespace std;

int main(int argc, char *argv[])
{
    GraphC h;
    h.read_gml(string(argv[1]));
//    h.read_weights_from_file(string(argv[2]));
    cout << h.is_weighted() << endl;

    CommunityStructure c(&h);
    c.set_random_seed();

    c.reindex_membership();
    c.sort_edges();
    try
    {
        c.optimize(QualityAsymptoticSurprise,MethodAgglomerative,1);
        c.print_membership();
        AsymptoticSurpriseFunction fun;
        cout << fun(h.get_igraph(),c.get_membership(),NULL) << endl;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }
    /*
    AgglomerativeOptimizer opt;
    opt.set_edges_order(c.get_sorted_edges_indices());

    Timer timer;timer.start();
    for (int i=0; i<100;++i)
        opt.optimize(h.get_igraph(),fun,c.get_membership());
    cout << timer.getElapsedTimeInMicroSec() << endl;
    c.reindex_membership();
    */
    //cout << "S=" << fun(h.get_igraph(),c.get_membership()) << endl;

    return 0;
}
