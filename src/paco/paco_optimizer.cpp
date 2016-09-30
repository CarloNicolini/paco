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

#include "Graph.h"
#include "Community.h"

using namespace std;



void exit_with_help()
{
//    std::printf(
//                "Usage: paco_optimizer [options] graph_file [output_suffix]\n"
//                "options:\n"
//                "-q [quality]"
//                "   0 Binary Surprise"
//                "   1 Significance"
//                "   2 Asymptotic Surprise"
//                "   3 Infomap\n"
//                "   4 Modularity\n"
//                "   5 Asymptotical Modularity (EXPERIMENTAL)\n"
//                "   6 Wonder (EXPERIMENTAL)\n"
//                "-m [method]:
//                "   0 Agglomerative Optimizer\n"
//                "   1 Random\n"
//                "   2 Simulated Annealing\n"
//                "-V [report_level] ERROR=0, WARNING=1, INFO=2, DEBUG=3, DEBUG1=4, DEBUG2=5, DEBUG3=6, DEBUG4=7\n"
//                "-S [seed] specify the random seed, default time(0)\n"
//                "-b [bool] wheter to start with initial random cluster or every node in its community\n"
//                "-r [repetitions], number of repetitions of PACO, default=1\n"
//                "-p [print solution]\n"
//                "\n"
//                );
    exit(1);
}



int main(int argc, char *argv[])
{
    srand(time(0));

    FILELog::ReportingLevel() = logDEBUG4;// static_cast<TLogLevel>(params.verbosity);

    //    Eigen::MatrixXd W(10,10);
    //    W << 0, 197, 0, 0, 101, 0, 0, 130, 0, 0, 197, 0, 0, 101, 109, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 138, 0, 164, 147, 0, 101, 0, 0, 0, 143, 151, 164, 167, 105, 101, 109, 0, 0, 0, 151, 159, 172, 130, 118, 0, 0, 130, 143, 151, 0, 0, 0, 0, 101, 0, 0, 138, 151, 159, 0, 0, 0, 0, 109, 130, 0, 0, 164, 172, 0, 0, 0, 0, 122, 0, 0, 164, 167, 130, 0, 0, 0, 0, 0, 0, 0, 147, 105, 118, 101, 109, 122, 0, 0;

    GraphC g;
    g.read_weighted_edge_list(string(argv[1]),51653);
    g.info();

    CommunityStructure comm(&g);
    comm.optimize(QualityAsymptoticSurprise,MethodAgglomerative);
    comm.reindex_membership();

    comm.save_membership("membership.csv",comm.get_membership());
    AsymptoticSurpriseFunction qual;
    cout << "Final quality=" << qual(g.get_igraph(),comm.get_membership(),g.get_edge_weights()) << endl;
    //    GraphC y;
    //    y.compute_vertex_degrees(false);
    //    y.compute_vertex_strenghts(false);

    //    GraphC x(y);

    return 0;
}

