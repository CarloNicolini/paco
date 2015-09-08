#include <iostream>
#include "AnnealOptimizer.h"
#include <set>
AnnealOptimizer::AnnealOptimizer(const igraph_t *g, const QualityFunction &fun,const  igraph_vector_t *memb, const igraph_vector_t *weights) : QualityOptimizer(g, fun, memb)
{
    this->optimize(g,fun,memb,weights);
}

AnnealOptimizer::~AnnealOptimizer()
{

}

double AnnealOptimizer::diff_move(const igraph_t *g, const QualityFunction &fun, const igraph_vector_t *memb, int vert, size_t dest_comm)
{
    // try to join the vertices
    double pre = fun(par);
    int orig_comm = memb->stor_begin[vert]; // save old original community of vert
    bool vertex_moved = par->move_vertex(g, memb,vert,dest_comm);
    if (vertex_moved)
    {
        double post = fun(par);
        if (post<pre)
        {
            par->move_vertex(g, memb,vert,orig_comm); // restore previous status
            return post-pre;
        }
        else
            return post-pre;
    }
    else
        return 0;
}

bool AnnealOptimizer::optimize(const igraph_t *g, const QualityFunction &fun, const  igraph_vector_t *memb, const igraph_vector_t *weights)
{
    par->init(g,memb);
    igraph_rng_seed(igraph_rng_default(), time(0));
    igraph_rng_t *rng = igraph_rng_default();
    int n = igraph_vcount(g);
    int nedges = igraph_ecount(g);
    size_t nhits = 0;
    double temp=this->param.temperature;
    size_t nstep = 0;
    size_t nImproves = 0;
    size_t nWorse = 0;

    //cerr << param.nIterations << " " << param.temperature << " " << param.temp_scale << " " << param.tolerance << " " << param.min_temp << endl;

    double best_val = 0;//std::numeric_limits<double>::min();
    igraph_vector_t best_memb;
    igraph_vector_init(&best_memb,n);
    // Start the optimization

    while (true)
    {
        ++nstep;
        double fval = fun(g,memb);
        if (fval > best_val)
        {
            best_val = fval;
            cout << best_val << endl;
            igraph_vector_update(&best_memb,memb);
        }

        temp = param.temperature*exp(-param.temp_scale*nstep/param.nIterations);
        // Choose a random edge
        int e = rand()%igraph_ecount(g);
        // Endpoints of random edge
        int ev1, ev2;
        igraph_edge(g,e,&ev1,&ev2);
        size_t cev2 = memb->stor_begin[ev2];
        // Get the difference in cost function of moving ev1 community to ev2 community
        double delta = diff_move(g,fun,memb,ev1,cev2);

        double rand_unif_01 = igraph_rng_get_unif01(rng);
        bool accept_better = (delta>0);
        bool accept_worse = exp(delta/temp) < rand_unif_01;
        if ( accept_better )
        {
            nImproves += size_t(accept_better);
        }
        else if ( accept_worse )
        {
            nWorse += size_t(accept_worse);
        }

        // Destroy with some random probability the current community, by selecting two random connected nodes and detaching them from the same community.
        if ( igraph_rng_get_unif01(rng) < 1E-3 )
        {
            int ec1=0,ec2=1;
            int count_comms = std::set<int>(memb->stor_begin,memb->stor_end).size();
            int c=0;
            while (c < count_comms)
            {
                igraph_edge(g,igraph_rng_get_integer(rng,0,nedges-1),&ec1,&ec2);
                if (memb->stor_begin[ec1] == memb->stor_begin[ec2])
                {
                    memb->stor_begin[ec1] = igraph_rng_get_integer(rng,0,n-1);
                }
                ++c;
            }
        }
        if ( fabs(delta)<param.tolerance)
            ++nhits;

        if (fval > param.minfval) // found a solution that is better than wanted solution
        {
            break;
        }
        else if (nhits > param.nHits || nstep>param.nIterations || temp < param.min_temp)
        {
            break;
        }
    }
    par->reindex(memb);
    // Copy the best solution to final membership
    for (size_t i=0; i<igraph_vector_size(memb); ++i)
        VECTOR(*memb)[i] = VECTOR(best_memb)[i];
    igraph_vector_destroy(&best_memb); // dispose memory
}
