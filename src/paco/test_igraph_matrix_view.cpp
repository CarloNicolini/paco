#include <iostream>
#include <string.h>
#include <sstream>
#include <igraph.h>
#include <sys/time.h>

#include "igraph_utils.h"
#include "Graph.h"
#include "Community.h"

#include "RandomOptimizer.h"
#include "AgglomerativeOptimizer.h"
#include "AnnealOptimizer.h"

#include "SurpriseFunction.h"
#include "AsymptoticSurpriseFunction.h"
#include "SignificanceFunction.h"

#include <Eigen/Core>

using namespace std;


void igraph_matrix_view2(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}


void create_igraph_t(const double *W, const int n, const int m, igraph_t *graph)
{
    igraph_matrix_t w_adj; // prepare the adjacency matrix
    if (n!=m)
    {
        throw std::logic_error("Non square adjacency matrix");
    }

    igraph_matrix_view(&w_adj,const_cast<igraph_real_t*>(W),n,m); // this aligns the igraph_matrix_t internal structures to the same addresses of `const double *W` matrix. No need to allocate new space.
    // Then create the igraph object
    igraph_weighted_adjacency(graph,&w_adj,IGRAPH_ADJ_UNDIRECTED,NULL,true);
}

int main()
{
    int n=100;
    Eigen::MatrixXd W;
    W.setRandom(n,n);
    W.array() += 1.0;

    for (int i=0; i<n; ++i)
        W(i,i)=0;

    W = W.triangularView<Eigen::Upper>();
    W = W.transpose().eval() + W;
    W(1,1)=1;
    //cout << W.topLeftCorner(10,10) << endl;
    GraphC *G;
    try
    {
        G = new GraphC(W);
        delete G;
    }
    catch (std::exception &e)
    {
        cerr << e.what() << endl;
    }

    return 0;

#ifdef OVERFLOW_TEST
    igraph_t G;
    //igraph_empty(&G,0,IGRAPH_UNDIRECTED);

    int n = 100;
    double *A;
    A = new double[n*n];
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
        {
            A[i*n+j] = A[j*n+i] = i+j;
            //A[j*n+i] = double(i+j);// to enforce symmetry
        }
    }

    // Set diagonal to zero
    for (int i=0; i<n; i++)
        A[i*n+i] = 0;

    igraph_matrix_t W;
    igraph_matrix_view(&W, A, n,n);

    create_igraph_t(A, n,n, &G);

    igraph_destroy(&G);
    delete[] A;
    return 0;
#endif
}
