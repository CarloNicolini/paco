#include <igraph/igraph.h>

int main()
{

    igraph_t g;
    const igraph_vector_t v=IGRAPH_VECTOR_NULL;
    igraph_real_t edges1[] = { 0,1, 1,2, 2,2, 2,3, 2,4, 3,4 };
    igraph_vector_t was;
    igraph_integer_t size;
    igraph_es_t it;
    long int i;

    igraph_vector_view(&v, edges1, sizeof(edges1)/sizeof(igraph_real_t));
    igraph_vector_init(&was, 0);

    /******************************************/
    /* Directed graph                         */
    /******************************************/

    igraph_create(&g, &v, 0, IGRAPH_DIRECTED);

    /* Simple test, all neighbors */
    for (i=0; i<=igraph_vector_max(&v); i++)
    {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_ALL);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it))
        {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    /* Simple test, outgoing neighbors */
    for (i=0; i<=igraph_vector_max(&v); i++)
    {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_OUT);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it))
        {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }


    /* Simple test, incoming neighbors */
    for (i=0; i<=igraph_vector_max(&v); i++)
    {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_IN);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it))
        {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    igraph_destroy(&g);

    /******************************************/
    /* Undirected graph                       */
    /******************************************/

    igraph_create(&g, &v, 0, IGRAPH_UNDIRECTED);

    /* Simple test, all neighbors */
    for (i=0; i<=igraph_vector_max(&v); i++)
    {
        igraph_vector_clear(&was);
        igraph_es_adj(&g, &it, i, IGRAPH_ALL);
        igraph_es_size(&g, &it, &size);
        printf("%ld\n", (long)size);
        while (!igraph_es_end(&g, &it))
        {
            igraph_vector_push_back(&was, igraph_es_adj_vertex(&g, &it));
            igraph_es_next(&g, &it);
        }
        igraph_es_destroy(&it);
        igraph_vector_sort(&was);
        igraph_vector_print(&was);
    }

    igraph_destroy(&g);
    igraph_vector_destroy(&was);

    return 0;
}
