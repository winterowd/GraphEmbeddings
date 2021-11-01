#include "AuxiliaryRoutinesForNauty.h"

/// set the partition for two-rooted graphs only1
void SetColoredPartition(int* c, int* lab, int* ptn, int n)
{
    /* loop over all colors to fill it */
    int cur_index = 0;
    for (int i = 0; i<3; i++) /// loop over colors
    {
        for (int j = 0; j<n; j++) /// loop over vertex labels
        {
            if (c[j] == i)
            {
                lab[cur_index] = j;
                ptn[cur_index] = 1;
                cur_index++;
            }
        }

        if (cur_index > 0) /// partition ends with a zero
        {
            ptn[cur_index - 1] = 0;
        }
    }
}

/** DERRICK STOLEE
 * SetColoredPartition sets up the nauty data structures lab & ptn according to a
 * given vertex coloring. After the execution, lab & ptn will be ready to send to nauty.
 *
 * @param n the number of vertices.
 * @param k the number of colors.
 * @param c the colors as in outputEColoredCanonicalLabeling.
 * @param lab an integer array of at least g->nv positions. WILL BE MODIFIED.
 * @param ptn an integer array of at least g->nv positions. WILL BE MODIFIED.
 */
void SetColoredPartition(int n, int k, int* c, int* lab, int* ptn)
{
    /* loop over all colos to fill it */
    int cur_index = 0;
    int i,j;
    for ( i = 0; i < k; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            if ( c[j] == i )
            {
                lab[cur_index] = j;
                ptn[cur_index] = 1;
                cur_index++;
            }
        }

        if (cur_index > 0)
        {
            ptn[cur_index - 1] = 0;
        }
    }
}

/**
 * SetColoredPartition sets up the nauty data structures lab & ptn according to a
 * given vertex coloring. After the execution, lab & ptn will be ready to send to nauty.
 *
 * @param g the graph as a sparsegraph.
 * @param k the number of colors.
 * @param c the colors as in outputEColoredCanonicalLabeling.
 * @param lab an integer array of at least g->nv positions. WILL BE MODIFIED.
 * @param ptn an integer array of at least g->nv positions. WILL BE MODIFIED.
 */
void SetColoredPartition(sparsegraph* g, int k, int* c, int* lab, int* ptn)
{
    /* loop over all colos to fill it */
    int cur_index = 0;
    int i,j;
    for ( i = 0; i < k; i++ )
    {
        for ( j = 0; j < g->nv; j++ )
        {
            if ( c[j] == i )
            {
                lab[cur_index] = j;
                ptn[cur_index] = 1;
                cur_index++;
            }
        }

        if (cur_index > 0)
        {
            ptn[cur_index - 1] = 0;
        }
    }
}

/// set up the colors based on the selected rooted vertices
/// i.e. rootedVertices[i] gets color i!
/// @param c: c-style int array of n elements
/// @param rootedVertices: array of vertices to be rooted
/// @param n: number of vertices
/// @param verbose: option for verbose output
void SetVertexColors(int *c, const std::vector<int>& rootedVertices, int n, bool verbose)
{
    if (rootedVertices.size()>n)
        throw std::invalid_argument("SetVertexColors requires number of rooted vertices to be less than or equal to N!\n");

    if (std::find_if(rootedVertices.begin(), rootedVertices.end(), [n](const int& x) { return x>=n; } ) != rootedVertices.end())
        throw std::invalid_argument("SetVertexColors requires each rooted vertex label x to satisfy 0 <= x < N!\n");

    if (verbose)
        std::cout << "SetVertexColors: " << rootedVertices.size() << " rooted vertices!\n";

    for (int i=0; i<n; ++i)
        c[i] = rootedVertices.size(); /// all vertices given color k

    for (int i=0; i<rootedVertices.size(); ++i) /// rooted vertices given colors 0,1,...,k-1
    {
        c[rootedVertices[i]] = i;
        if (verbose)
            std::cout << "SetVertexColors: vertex " << rootedVertices[i] << " assigned color " << i << "\n";
    }
}
