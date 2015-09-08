#include "Graph.h"

using namespace std;

#include <Eigen/Core>

int main(int argc, char *argv[])
{
    srand(time(0));

    FILELog::ReportingLevel() = logDEBUG4;// static_cast<TLogLevel>(params.verbosity);

    Eigen::MatrixXd W(10,10);
    W << 0, 197, 0, 0, 101, 0, 0, 130, 0, 0, 197, 0, 0, 101, 109, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 138, 0, 164, 147, 0, 101, 0, 0, 0, 143, 151, 164, 167, 105, 101, 109, 0, 0, 0, 151, 159, 172, 130, 118, 0, 0, 130, 143, 151, 0, 0, 0, 0, 101, 0, 0, 138, 151, 159, 0, 0, 0, 0, 109, 130, 0, 0, 164, 172, 0, 0, 0, 0, 122, 0, 0, 164, 167, 130, 0, 0, 0, 0, 0, 0, 0, 147, 105, 118, 101, 109, 122, 0, 0;

    GraphC *g = new GraphC();
    g->read_gml(string(argv[1]));
    g->read_gml(string(argv[1]));
    g->read_gml(string(argv[1]));
    g->read_weights_from_file(string(argv[2]));
    g->info();
    delete g;

    return 0;
}
