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


#include "Graph.h"
#include "Community.h"

using namespace std;



void exit_with_help()
{
    std::printf(
                "Usage: paco_optimizer graph_file [options]\n"
                "graph_file the file containing the graph. Accepted formats are pajek, graph_ml, adjacency matrix or"
                "\nedges list (the ncol format), additionally with a third column with edge weights"
                "options:\n"
                "-q [quality]\n"
                "   0 Binary Surprise\n"
                "   1 Significance\n"
                "   2 Asymptotic Surprise\n"
                "   3 Infomap\n"
            #ifdef EXPERIMENTAL_FEATURES
                "   4 Modularity (EXPERIMENTAL)\n"
                "   5 Asymptotical Modularity (EXPERIMENTAL)\n"
                "   6 Wonder (EXPERIMENTAL)\n"
            #endif
                "-m [method]:"
                "   0 Agglomerative Optimizer\n"
                "   1 Random\n"
                "   2 Simulated Annealing\n"
                "-V [report_level] ERROR=0, WARNING=1, INFO=2, DEBUG=3, DEBUG1=4, DEBUG2=5, DEBUG3=6, DEBUG4=7\n"
                "-S [seed] specify the random seed, default time(0)\n"
                "-b [bool] wheter to start with initial random cluster or every node in its community\n"
                "-r [repetitions], number of repetitions of PACO, default=1\n"
                "-p [print solution]\n"
                "\n"
                );
    exit(1);
}

enum error_type
{
    NO_ERROR = 0,
    ERROR_TOO_MANY_OUTPUT_ARGS = 1,
    ERROR_NOT_ENOUGH_ARGS = 2,
    ERROR_ARG_VALUE = 3,
    ERROR_ARG_TYPE = 4,
    ERROR_MATRIX = 5,
    ERROR_ARG_EMPTY=6,
    ERROR_ARG_UNKNOWN=7
};

static const char *error_strings[] =
{
    "",
    "Too many output arguments.",
    "Not enough input arguments.",
    "Non valid argument value.",
    "Non valid argument type.",
    "Non valid input adjacency matrix. PACO accepts symmetric real dense-type (n x n) matrices or sparse edges-list representation \
    [num_edges x 3] array of edges list with edge endpoints and weight.",
    "Expected some argument value but empty found.",
    "Unkwown argument."
};

struct PacoParams
{
    OptimizerType method=MethodAgglomerative;
    QualityType qual=QualitySurprise;
    size_t nrep=1;      // Maximum number of consecutive repetitions to perform.
    int rand_seed=-1; // random seed for the louvain algorithm
    int verbosity_level=0;
    std::string membership_file="membership.txt";
    std::string filename="";
    bool print_info=false;
};

/**
 * @brief parse_command_line
 * @param argc
 * @param argv
 * @param input_file_name
 */
PacoParams parse_command_line(int argc, char **argv)
{
    PacoParams params;

    int i=1;
    for(i=1; i<argc; i++)
    {
        if(argv[i][0] != '-')
            break;
        if(++i>=argc)
            exit_with_help();
        switch(argv[i-1][1])
        {
        case 'v':
        case 'V':
        {
            params.verbosity_level = atoi(argv[i]);
            if (params.verbosity_level>7)
                params.verbosity_level=7;
            FILELog::ReportingLevel() =  static_cast<TLogLevel>(params.verbosity_level);
            break;
        }
        case 's':
        case 'S':
        {
            params.rand_seed = atoi(argv[i]);
            break;
        }
        case 'q':
        case 'Q':
        {
            switch (atoi(argv[i]))
            {
            case 0:
            {

                params.qual = QualitySurprise;
                break;
            }
            case 1:
            {
                params.qual = QualitySignificance;
                break;
            }
            case 2:
            {
                params.qual = QualityAsymptoticSurprise;
                break;
            }
            case 3:
            {
                params.qual = QualityInfoMap;
                break;
            }
            default:
            {
                exit_with_help();
            }
            }
            break;
        }
        case 'm':
        case 'M':
        {
            params.method = static_cast<OptimizerType>(atoi(argv[i]));
            break;
        }
        case 'r':
        case 'R':
        {
            params.nrep = atoi(argv[i]);
            break;
        }
        case 'o':
        case 'O':
        {
            params.membership_file = std::string(argv[i]);
            break;
        }
        case 'p':
        case 'P':
        {
            params.print_info = (bool)atoi(argv[i]);
            break;
        }
        default:
            FILE_LOG(logERROR) << "Unknown option: " << argv[i-1][1] ;
            exit_with_help();
        }
    }

    // Determine filenames
    if(i>=argc)
        exit_with_help();

    char input[1024];
    strcpy(input, argv[i]);
    params.filename = std::string(input);

    std::ifstream is(input);
    if (!is.good())
    {
        cout << std::string("File \"" + params.filename + "\" not found") << endl;
        exit_with_help();
    }
    return params;
}


int main(int argc, char *argv[])
{
    PacoParams pars = parse_command_line(argc,argv);

    // Determine file format
    GraphC g;
    if(pars.filename.substr(pars.filename.find_last_of(".") + 1) == "net")
        g.read_pajek(pars.filename);
    else if(pars.filename.substr(pars.filename.find_last_of(".") + 1) == "adj")
        g.read_adj_matrix(pars.filename);
    else if(pars.filename.substr(pars.filename.find_last_of(".") + 1) == "gml")
        g.read_gml(pars.filename);
    else if (pars.filename.substr(pars.filename.find_last_of(".") + 1) == "ncol" || pars.filename.substr(pars.filename.find_last_of(".") + 1) == "edge")
        g.read_edge_list(pars.filename);
    else if(pars.filename.substr(pars.filename.find_last_of(".") + 1) == "wncol")
        g.read_weighted_edge_list(pars.filename);
    else
    {
        cerr << "Non supported graph format" << endl;
        exit_with_help();
    }

    if (pars.print_info)
    {
        FILELog::ReportingLevel() =  TLogLevel::logINFO;
        g.info();
    }

    CommunityStructure comm(&g);
    comm.set_random_seed(pars.rand_seed);
    double quality = comm.optimize(pars.qual,pars.method,pars.nrep);
    comm.reindex_membership();
    comm.save_membership(pars.membership_file.c_str(),comm.get_membership());
    cout << quality << endl;

    return 0;
}
