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
#include <fstream>

#include <igraph.h>
#include <igraph_iterators.h>

#include "Graph.h"
#include "SurpriseFunction.h"
#include "SignificanceFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "RandomOptimizer.h"
#include "AnnealOptimizer.h"
#include "Community.h"

#include "PartitionHelper.h"
#include "AgglomerativeOptimizer.h"

using namespace std;

int main(int argc, char *argv[])
{

    srand(time(0));

    FILELog::ReportingLevel() = logDEBUG4;// static_cast<TLogLevel>(params.verbosity);
    GraphC g;
    //g.read_gml("/home/carlo/workspace/BCT/Coactivation_matrix.gml");
    g.read_gml(argv[1]);

    CommunityStructure c(&g);
    c.read_membership_from_file(argv[2]);


    QualityFunction *fun;
    fun = dynamic_cast<SurpriseFunction*>( new SurpriseFunction);
    cout << (*fun)(g.get_igraph(),c.get_membership()) << endl;
    delete fun;

    QualityOptimizer *opt;
    opt = ( new AgglomerativeOptimizer);

    delete fun;
    delete opt;

    return 0;
}
