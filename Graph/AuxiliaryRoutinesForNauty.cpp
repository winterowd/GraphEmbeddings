#include "AuxiliaryRoutinesForNauty.h"

/// set the partition for two-rooted graphs only
void AuxiliaryRoutinesForNauty::SetColoredPartition(int* c, int* lab, int* ptn, int n)
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
void AuxiliaryRoutinesForNauty::SetColoredPartition(int n, int k, int* c, int* lab, int* ptn)
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
void AuxiliaryRoutinesForNauty::SetColoredPartition(sparsegraph* g, int k, int* c, int* lab, int* ptn)
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
void AuxiliaryRoutinesForNauty::SetVertexColors(int *c, const std::vector<int>& rootedVertices, int n, bool verbose)
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

/// routine which takes a graph in dense nauty format and returns the canonical graph
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(int n, graph *g)
{
    int m = SETWORDSNEEDED(n);
    std::string g6String = ntog6(g,m,n); /// graph to g6 string routine
    g6String.erase(std::remove(g6String.begin(), g6String.end(), '\n'), g6String.end()); // get rid of newline character
    return GetCanonicalGraphNauty(n, g6String);
}

/// routine which takes a g6 string and rooted vertices and returns the canonical graph (useful for debugging purposes etc.)
/// @arg n: number of vertices
/// @arg g6string: valid g6 string (user must check this!)
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(int n, const std::string& g6String)
{
    char *tempg6 = new char[g6String.length()+1];
    std::strcpy(tempg6, g6String.c_str());

    if (n != graphsize(tempg6)) /// check size
        throw std::invalid_argument("GetCanonicalGraph g6String not of size N!\n");

    int mWords = SETWORDSNEEDED(n); /// set m

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    DYNALLOC2(graph, g, g_sz, n, mWords, "malloc"); /// allocate graph
    DYNALLOC2(graph, cg, cg_sz, n, mWords, "malloc"); /// allocate canonical graph
    DYNALLOC2(int, lab, lab_sz, n, mWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, mWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, mWords, "malloc");

    stringtograph(tempg6, g, mWords); /// g6 string to densenauty
    GraphContainer resultContainer(n, mWords, g); /// container for graphs

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.getcanon = true; /// get canong

    /// call densenauty
    densenauty(g, lab, ptn, orbits, &options, &stats, mWords, n, cg);

    resultContainer.CanonicalRelabeling(lab);

#ifdef DEBUG
    GraphContainer refContainer(n, mWords, cg);
    if (refContainer==refContainer)
        std::cout << "Success! resultContainer equals refContainer!\n";
    else
        std::cout << "ERROR: resultContainer does not equal refContainer!\n";
#endif
    delete[] tempg6;

    DYNFREE(g, g_sz); /// free graph
    DYNFREE(cg, cg_sz); /// free canonical graph
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    return resultContainer;
}

/// routine which takes a g6 string and rooted vertices and returns the canonical colored graph (useful for debugging purposes etc.)
/// @arg n: number of vertices
/// @arg g6string: valid g6 string (user must check this!)
/// @arg rootedVertices: labels of rooted (colored) vertices (BEFORE relabeling!)
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(int n, const std::string& g6String, const std::vector<int>& rootedVertices)
{
    char *tempg6 = new char[g6String.length()+1];
    std::strcpy(tempg6, g6String.c_str());

    if (n != graphsize(tempg6)) /// check size
        throw std::invalid_argument("GetCanonicalColoredGraph g6String not of size N!\n");

    if (rootedVertices.size() != 2)
        throw std::invalid_argument("GetCanonicalColoredGraph rootedVertices should be of size 2!\n");

    if (rootedVertices[0] < 0 || rootedVertices[0] >= n)
        throw std::invalid_argument("GetCanonicalColoredGraph requires 0 <= rootedVertices[0] < N!\n");

    if (rootedVertices[1] < 0 || rootedVertices[1] >= n)
        throw std::invalid_argument("GetCanonicalColoredGraph requires 0 <= rootedVertices[1] < N!\n");

    int mWords = SETWORDSNEEDED(n); /// set m

    int *c = (int*)malloc(n * sizeof(int)); /// alloc C-style arrray for colors of vertices

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    DYNALLOC2(graph, g, g_sz, n, mWords, "malloc"); /// allocate graph
    DYNALLOC2(graph, cg, cg_sz, n, mWords, "malloc"); /// allocate canonical graph
    DYNALLOC2(int, lab, lab_sz, n, mWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, mWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, mWords, "malloc");

    stringtograph(tempg6, g, mWords); /// g6 string to densenauty
    GraphContainer resultContainer(n, mWords, g, true, 2); /// container for graphs

    resultContainer.SetRootedVertex(0, rootedVertices[0]);
    resultContainer.SetRootedVertex(1, rootedVertices[1]);

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canong

    SetVertexColors(c, rootedVertices, n);
    SetColoredPartition(c, lab, ptn, n); /// set coloring (will be written over in call to densenauty)

    /// call densenauty
    densenauty(g, lab, ptn, orbits, &options, &stats, mWords, n, cg);

    resultContainer.ColoredCanonicalRelabeling(lab, rootedVertices[0], rootedVertices[1]);

#ifdef DEBUG
    GraphContainer refContainer(n, mWords, cg, true, 2);
    refContainer.SetRootedVertex(0, 0);
    refContainer.SetRootedVertex(1, 1);
    if (refContainer==refContainer)
        std::cout << "Success! resultContainer equals refContainer!\n";
    else
        std::cout << "ERROR: resultContainer does not equal refContainer!\n";
#endif

    free(c); /// free vertex colors

    delete[] tempg6;

    DYNFREE(g, g_sz); /// free graph
    DYNFREE(cg, cg_sz); /// free canonical graph
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    return resultContainer;
}

