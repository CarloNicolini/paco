#include "benchm.h"
#include "FileLogger.h"

/**
 * @brief main
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char * argv[])
{
    FILELog::ReportingLevel() = logDEBUG;
    srand_file();
    Parameters p;
    if(set_parameters(argc, argv, p)==false)
    {
        if (argc>1)
            cerr << "Please, look at README.txt..." << endl;
        return -1;
    }

    erase_file_if_exists("network.dat");
    erase_file_if_exists("matrix.dat");
    erase_file_if_exists("community.dat");
    erase_file_if_exists("statistics.dat");

    Eigen::MatrixXd W; // weighted adjacency matrix
    vector<int> m(p.num_nodes);
    benchmark(p.excess, p.defect, p.num_nodes, p.average_k, p.max_degree, p.tau, p.tau2, p.mixing_parameter,  p.mixing_parameter2,  p.beta, p.overlapping_nodes, p.overlap_membership, p.nmin, p.nmax, p.fixed_range, p.clustering_coeff, W,m);

    ofstream os; os.open("matrix.dat");
    os << W << endl;

    return 0;
}
