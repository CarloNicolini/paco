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

using namespace std;

int main(int argc, char *argv[])
{
    GraphC h;
    h.read_adj_matrix(string(argv[1]));
    const igraph_t *g = h.get_igraph();
    CommunityStructure c(&h);
    c.read_membership_from_file(argv[2]);

    const igraph_vector_t *m = c.get_membership();
    const igraph_vector_t *w = h.get_edge_weights();
    PartitionHelper par;
    par.init(g,m,w);

    AsymptoticSurpriseFunction f;
    cout << "AS pre movement=" << f(g,m,w) << endl;
    par.move_vertex(g,m,0,29,h.get_edge_weights());
    par.move_vertex(g,m,1,29,h.get_edge_weights());
    par.move_vertex(g,m,2,29,h.get_edge_weights());
    par.move_vertex(g,m,3,29,h.get_edge_weights());
    //igraph_vector_print(m);
    par.print();
    cout << "AS pre movement=" << f(g,m,w) << endl;
    /*
    SurpriseFunction f;
    cout << f(&par) << endl;

    for (int i=0; i<6; ++i)
    {
        par.move_vertex(g,m,i,11);
        cout << f(&par) << endl;
    }

    for (int i=0; i<6; ++i)
    {
        par.move_vertex(g,m,i,6);
        cout << f(&par) << endl;
    }
*/

    return 0;
}
