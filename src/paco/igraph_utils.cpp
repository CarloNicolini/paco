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


#include <algorithm>
#include "igraph_utils.h"

void igraph_matrix_view(igraph_matrix_t *A, igraph_real_t *data, int nrows, int ncols)
{
    A->ncol = nrows;
    A->nrow = ncols;
    A->data.stor_begin = data;
    A->data.stor_end = data+ncols*nrows;
    A->data.end = A->data.stor_end;
}

igraph_vector_t* order_membership(const igraph_vector_t *curmemb)
{
    int nVertices = igraph_vector_size(curmemb);
    // groupmap
    igraph_vector_t group_map; igraph_vector_init(&group_map,nVertices);
    igraph_vector_fill(&group_map,-1);
    int current_number=0;
    for (int i=0; i<nVertices; ++i)
    {
        if (group_map.stor_begin[int(curmemb->stor_begin[i])] < 0)
        {
            group_map.stor_begin[int(curmemb->stor_begin[i])] = current_number;
            ++current_number;
        }
    }


    igraph_vector_t *newmemb = new igraph_vector_t;
    igraph_vector_init(newmemb,nVertices);
    igraph_vector_fill(newmemb,-1);

    for (int i=0; i<nVertices; ++i)
    {
        int val = group_map.stor_begin[int(curmemb->stor_begin[i])];
        if (val<0)
            val = *(group_map.stor_end);
        newmemb->stor_begin[i] = val;
    }


    // Clean memory
    igraph_vector_destroy(&group_map);

    return newmemb;
}

// Implementation based on igraph original functions of the index of structural similarity described in
// "Density-based shrinkage for revealing hierarchical and overlapping community structure in networks"
int igraph_i_neisets_intersect(const igraph_t *graph, const igraph_vector_t *v1, const igraph_vector_t *v2, const igraph_vector_t *weights, double *weight_union, double *weight_intersection)
{
    std::vector<int> vv1 = std::vector<int>(v1->stor_begin,v1->stor_end);
    std::vector<int> vv2 = std::vector<int>(v2->stor_begin,v2->stor_end);

    std::vector<int> vunion,vintersection;
    std::set_intersection(vv1.begin(),vv1.end(),vv2.begin(),vv2.end(),  std::back_inserter(vintersection));
    std::set_union(vv1.begin(),vv1.end(),vv2.begin(),vv2.end(),         std::back_inserter(vunion));

    /* ASSERT: v1 and v2 are sorted */
    //long int i=0, j=0, v1size=igraph_vector_size(v1), v2size=igraph_vector_size(v2);
    *weight_union = 0;
    *weight_intersection = 0;
    /*
    while (i < v1size && j < v2size)
    {
        if (VECTOR(*v1)[i] == VECTOR(*v2)[j])
        {
            int x=VECTOR(*v1)[i];
            int y=VECTOR(*v2)[j];
            int e=igraph_get_eid(graph,&e,x,y,false,false);
            (*weight_intersection) += VECTOR(*weights)[e];
            (*weight_union)-=VECTOR(*weights)[e];
            i++;
            j++;
        }
        else if (VECTOR(*v1)[i] < VECTOR(*v2)[j])
            i++;
        else
            j++;
    }
*/
    //*weight_intersection=0;
    //*weight_intersection=0;
    return 0;
}

int igraph_similarity_jaccard_weighted_pairs(const igraph_t *graph, igraph_vector_t *res, const igraph_vector_t *pairs, const igraph_vector_t *weights, igraph_neimode_t mode, igraph_bool_t loops)
{
    igraph_lazy_adjlist_t al;
    long int i;
    long int j;
    long int k;
    long int u;
    long int v;

    igraph_vector_t *v1, *v2;

    k = igraph_vector_size(pairs);
    if (k % 2 != 0)
        IGRAPH_ERROR("number of elements in `pairs' must be even", IGRAPH_EINVAL);
    IGRAPH_CHECK(igraph_vector_resize(res, k/2));

    IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &al, mode, IGRAPH_SIMPLIFY));
    IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &al);
    /*
    igraph_bool_t *seen;
    if (loops)
    {
        // Add the loop edges
        i = igraph_vcount(graph);
        seen = igraph_Calloc(i, igraph_bool_t);
        if (seen == 0)
            IGRAPH_ERROR("cannot calculate Jaccard similarity", IGRAPH_ENOMEM);
        IGRAPH_FINALLY(free, seen);

        for (i = 0; i < k; i++) {
            j = (long int) VECTOR(*pairs)[i];
            if (seen[j])
                continue;
            seen[j] = 1;
            v1=igraph_lazy_adjlist_get(&al, (igraph_integer_t) j);
            if (!igraph_vector_binsearch(v1, j, &u))
                igraph_vector_insert(v1, u, j);
        }

        free(seen);
        IGRAPH_FINALLY_CLEAN(1);
    }
*/
    for (i = 0, j = 0; i < k; i += 2, j++)
    {
        u = (long int) VECTOR(*pairs)[i];
        v = (long int) VECTOR(*pairs)[i+1];

        if (u == v)
        {
            VECTOR(*res)[j] = 1.0;
            continue;
        }

        v1 = igraph_lazy_adjlist_get(&al, (igraph_integer_t) u);
        v2 = igraph_lazy_adjlist_get(&al, (igraph_integer_t) v);

        /********************/
        std::vector<int> vv1 = std::vector<int>(v1->stor_begin,v1->stor_end);
        std::vector<int> vv2 = std::vector<int>(v2->stor_begin,v2->stor_end);
        vv1.push_back(u);
        vv2.push_back(v);

        std::vector<int> vunion,vintersection;
        std::set_intersection(vv1.begin(),vv1.end(),vv2.begin(),vv2.end(),  std::back_inserter(vintersection));
        std::set_union(vv1.begin(),vv1.end(),vv2.begin(),vv2.end(),         std::back_inserter(vunion));


        double weight_intersection=0;

        for (std::vector<int>::iterator x = vintersection.begin(); x!=vintersection.end();++x)
        {
            int e1=-1,e2=-1;
            igraph_get_eid(graph,&e1,u,*x,false,false);
            igraph_get_eid(graph,&e2,v,*x,false,false);
            double wux = weights->stor_begin[e1];
            double wvx = weights->stor_begin[e2];
            weight_intersection += wux*wvx;
        }

        double sumwux2=0;
        double sumwvx2=0;
        for (std::vector<int>::iterator x = vv1.begin(); x!=vv1.end(); ++x)
        {
            int e=-1;
            igraph_get_eid(graph,&e,u,*x,false,false);
            double wux2 = SQR(weights->stor_begin[e]);// * weights->stor_begin[e];
            sumwux2+=wux2;
        }

        for (std::vector<int>::iterator x = vv2.begin(); x!=vv2.end(); ++x)
        {
            int e=-1;
            igraph_get_eid(graph,&e,v,*x,false,false);
            double wvx2=SQR(weights->stor_begin[e]);
            sumwvx2+=wvx2;
        }

        double weight_union = sqrt(sumwux2)*sqrt(sumwvx2);

        VECTOR(*res)[j] = weight_intersection/weight_union;
    }

    igraph_lazy_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

int igraph_similarity_jaccard_weighted_es(const igraph_t *graph, igraph_vector_t *res,
                                          const igraph_es_t es, const igraph_vector_t *weights, igraph_neimode_t mode, igraph_bool_t loops)
{
    igraph_vector_t v;
    igraph_eit_t eit;

    IGRAPH_VECTOR_INIT_FINALLY(&v, 0);

    IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
    IGRAPH_FINALLY(igraph_eit_destroy, &eit);

    while (!IGRAPH_EIT_END(eit)) {
        long int eid = IGRAPH_EIT_GET(eit);
        igraph_vector_push_back(&v, IGRAPH_FROM(graph, eid));
        igraph_vector_push_back(&v, IGRAPH_TO(graph, eid));
        IGRAPH_EIT_NEXT(eit);
    }

    igraph_eit_destroy(&eit);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_similarity_jaccard_weighted_pairs(graph, res, &v, weights, mode, loops));
    igraph_vector_destroy(&v);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

int igraph_read_graph_weighted_edgelist(igraph_t *graph, FILE *instream,
                                        igraph_integer_t n, igraph_bool_t directed, igraph_vector_t *edge_weights)
{

    igraph_vector_t edges=IGRAPH_VECTOR_NULL;
    long int from, to;
    double weight;
    int c;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, 100));

    /* skip all whitespace */
    do {
        c = getc (instream);
    } while (isspace (c));
    ungetc (c, instream);

    int linenumber=1;
    std::vector<double> eweights;
    while (!feof(instream))
    {
        int read;

        //IGRAPH_ALLOW_INTERRUPTION();

        read=fscanf(instream, "%li", &from);
        if (read != 1)
        {
            IGRAPH_ERROR("parsing edgelist file failed", IGRAPH_PARSEERROR);
        }
        read=fscanf(instream, "%li", &to);
        if (read != 1)
        {
            IGRAPH_ERROR("parsing edgelist file failed", IGRAPH_PARSEERROR);
        }
        read=fscanf(instream, "%lf", &weight);
        if (read != 1)
        {
            IGRAPH_ERROR("parsing edgelist file failed (invalid weight)", IGRAPH_PARSEERROR);
        }
        IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
        IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
        eweights.push_back(weight);
        /* skip all whitespace */
        do {
            c = getc (instream);
        } while (isspace (c));
        ungetc (c, instream);
        ++linenumber;
    }
    // copy all the weights to the weight vector
    igraph_vector_init_copy(edge_weights, eweights.data(), eweights.size());

    IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * @brief free_complist See the example in documentation about igraph_decompose
 * @param complist
 */
void free_complist(igraph_vector_ptr_t *complist)
{
    long int i;
    for (i=0; i<igraph_vector_ptr_size(complist); i++)
    {
        igraph_t *c = static_cast<igraph_t*>(complist->stor_begin[i]);
        igraph_destroy(c);
        free(c);
    }
}
