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

using namespace std;

int main(int argc, char *argv[])
{
    GraphC h;
    h.read_gml(argv[1]);
    //h.read_weights_from_file(argv[2]);

    CommunityStructure c(&h);
    c.set_random_seed();

    //if (argc<=2)
      //  c.read_membership_from_file(argv[3]);
    c.reindex_membership();
    c.sort_edges();

    AgglomerativeOptimizer opt;
    AsymptoticSurpriseFunction fun;
    opt.set_edges_order(c.get_sorted_edges_indices());
    opt.optimize(h.get_igraph(),fun,c.get_membership());
    c.reindex_membership();

    cout << "S=" << fun(h.get_igraph(),c.get_membership()) << endl;

    return 0;
}
