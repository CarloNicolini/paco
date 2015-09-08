#include <iostream>

#include "Graph.h"
#include "SurpriseFunction.h"
#include "Community.h"
#include "PartitionHelper.h"

using namespace std;

int main(int argc, char *argv[])
{
    GraphC h;
    h.read_gml(argv[1]);
    const igraph_t *g = h.get_igraph();
    CommunityStructure c(&h);
    c.read_membership_from_file("../data/test.csv");
    c.reindex_membership();

    const igraph_vector_t *m = c.get_membership();
    PartitionHelper par;
    par.init(g,m);

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

    return 0;
}
