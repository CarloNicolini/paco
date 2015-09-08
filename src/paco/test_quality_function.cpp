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

    return 0;
}
