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

#include <Eigen/Core>
#include "Graph.h"
#include "Community.h"

using namespace std;

vector<double> read_vector(const string &filename)
{
    std::ifstream membership_file(filename.c_str());
    if (!membership_file.good())
        throw std::runtime_error(std::string("File not found"));

    string line;
    vector<igraph_real_t> vmemb;
    while ( getline(membership_file,line) )
    {
        if (line.size() == 0)
            continue;
        std::stringstream str(line);
        igraph_real_t memb;
        str >> memb;
        vmemb.push_back(memb);
    }
    return vmemb;
}


int main(int argc, char *argv[])
{
    srand(time(0));
    FILELog::ReportingLevel() = logDEBUG4;// static_cast<TLogLevel>(params.verbosity);
    
    vector<double> i = read_vector("i.csv");
    vector<double> j = read_vector("j.csv");
    vector<double> w = read_vector("w.csv");

    vector<double> e;
    for (int l=0; l<i.size(); ++l)
    {
        e.push_back(i[l]-1);
        e.push_back(j[l]-1);
    }

    GraphC *g = new GraphC(e.data(),w.data(),w.size());

    g->info();
    g->print();

//    CommunityStructure c(g);
//    //c.set_random_seed(pars.rand_seed);
//    double finalquality=c.optimize(QualitySurprise,MethodAgglomerative,5);
//    cout << finalquality << endl;
    delete g;
    return 0;
}
