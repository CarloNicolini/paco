#include "benchm.h"

int main(int argc, char * argv[])
{
    srand_file();
    Parameters p;
    if(set_parameters(argc, argv, p)==false)
    {
        if (argc>1)
            cerr<<"Please, look at ReadMe.txt..."<<endl;
        return -1;
    }

    erase_file_if_exists("network.dat");
    erase_file_if_exists("community.dat");
    erase_file_if_exists("statistics.dat");

    Eigen::MatrixXd W; // weighted adjacency matrix
    vector<int> m(p.num_nodes);
    benchmark(p.excess, p.defect, p.num_nodes, p.average_k, p.max_degree, p.tau, p.tau2, p.mixing_parameter,  p.mixing_parameter2,  p.beta, p.overlapping_nodes, p.overlap_membership, p.nmin, p.nmax, p.fixed_range, p.clustering_coeff, W,m);
    //cerr << W.block(0,0,20,20) << endl;
    return 0;
}
