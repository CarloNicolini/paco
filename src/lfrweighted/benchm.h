#include "binary_benchm.h"

int print_network(deque<set<int> > & E, const deque<deque<int> > & member_list, const deque<deque<int> > & member_matrix,
                  deque<int> & num_seq, deque<map <int, double > > & neigh_weigh, double beta, double mu, double mu0);
int check_weights(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list,
                  deque<deque<double> > & wished, deque<deque<double> > & factual, const double tot_var, double *strs);
int propagate(deque<map <int, double > > & neigh_weigh, const deque<deque<int> > & member_list,
              deque<deque<double> > & wished, deque<deque<double> > & factual, int i, double & tot_var, double *strs, const deque<int> & internal_kin_top);
int weights(deque<set<int> > & en, const deque<deque<int> > & member_list, const double beta, const double mu, deque<map <int, double > > & neigh_weigh);
int benchmark(bool excess, bool defect, int num_nodes, double  average_k, int  max_degree, double  tau, double  tau2, double  mixing_parameter, double  mixing_parameter2, double  beta, int  overlapping_nodes, int  overlap_membership, int  nmin, int  nmax, bool  fixed_range, double ca);
void erase_file_if_exists(string s);
