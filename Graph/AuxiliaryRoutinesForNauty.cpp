#include "AuxiliaryRoutinesForNauty.h"

/// set the partition for one- and two-rooted graphs only
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
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(int n, graph *g, int* labOutput)
{
    int m = SETWORDSNEEDED(n);
    std::string g6String = ntog6(g,m,n); /// graph to g6 string routine
    g6String.erase(std::remove(g6String.begin(), g6String.end(), '\n'), g6String.end()); // get rid of newline character
    return GetCanonicalGraphNauty(n, g6String, labOutput);
}

/// routine which takes a g6 string and rooted vertices and returns the canonical graph (useful for debugging purposes etc.)
/// @arg n: number of vertices
/// @arg g6string: valid g6 string (user must check this!)
/// @arg labOutput: labeling for output (if NULL, allocate locally; assume it is properly allocated by calling routine otherwise)
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(int n, const std::string& g6String, int *labOutput)
{
    char *tempg6 = new char[g6String.length()+1];
    std::strcpy(tempg6, g6String.c_str());

    if (n != graphsize(tempg6)) /// check size
        throw std::invalid_argument("GetCanonicalGraphNauty g6String not of size N!\n");

    int mWords = SETWORDSNEEDED(n); /// set m
    bool allocateLabLocally = (labOutput==NULL) ? true : false; /// do we need to allocate memory for label?

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    DYNALLOC2(graph, g, g_sz, n, mWords, "malloc"); /// allocate graph
    DYNALLOC2(graph, cg, cg_sz, n, mWords, "malloc"); /// allocate canonical graph
    size_t lab_sz=0;
    if (allocateLabLocally)
        DYNALLOC2(int, labOutput, lab_sz, n, mWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, mWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, mWords, "malloc");

    stringtograph(tempg6, g, mWords); /// g6 string to densenauty
    GraphContainer resultContainer(n, mWords, g); /// container for graphs

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.getcanon = true; /// get canong

    /// call densenauty
    densenauty(g, labOutput, ptn, orbits, &options, &stats, mWords, n, cg);

    if (stats.errstatus!=0)
        std::cerr << "AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty: call to NAUTY produced an error!\n";

    resultContainer.CanonicalRelabeling(labOutput); /// relabel

    resultContainer.SetSymmFactor(AuxiliaryRoutinesForNauty::ExtractSymmFactor(stats)); /// get symmetry factor

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
    if (allocateLabLocally) /// deallocate memory and set pointer back to NULL
    {
        DYNFREE(labOutput,lab_sz);
        labOutput = NULL;
    }
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    return resultContainer;
}

/// routine which takes a g6 string and rooted vertices and returns the canonical colored graph (useful for debugging purposes etc.)
/// @arg n: number of vertices
/// @arg g6string: valid g6 string (user must check this!)
/// @arg rootedVertices: labels of rooted (colored) vertices (BEFORE relabeling!) (NOTE: only one or two rooted vertices allowed!) (labeled 0, 1, ..., N-1 !!!!!)
/// @arg labOutput: labeling for output (if NULL, allocate locally; assume it is properly allocated by calling routine otherwise)
GraphContainer AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(int n, const std::string& g6String, const std::vector<int>& rootedVertices, int *labOutput)
{
    char *tempg6 = new char[g6String.length()+1];
    std::strcpy(tempg6, g6String.c_str());

    if (n != graphsize(tempg6)) /// check number of vertices
        throw std::invalid_argument("GetCanonicalColoredGraphNauty g6String not of size N!\n");

    if (rootedVertices.size() != 1 && rootedVertices.size() != 2) /// check number of rooted vertices
        throw std::invalid_argument("GetCanonicalColoredGraphNauty rootedVertices should be of size 1 or 2!\n");

    for (int i=0; i<rootedVertices.size(); ++i) /// check labels of rooted vertices
        if (rootedVertices[i] < 0 || rootedVertices[i] >= n)
            throw std::invalid_argument("GetCanonicalColoredGraphNauty requires 0 <= rootedVertices[i] < N!\n");

    bool allocateLabLocally = (labOutput==NULL) ? true : false; /// do we need to allocate memory for label?
    int mWords = SETWORDSNEEDED(n); /// set m

    int *c = (int*)malloc(n * sizeof(int)); /// alloc C-style arrray for colors of vertices

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    DYNALLOC2(graph, g, g_sz, n, mWords, "malloc"); /// allocate graph
    DYNALLOC2(graph, cg, cg_sz, n, mWords, "malloc"); /// allocate canonical graph
    size_t lab_sz = 0;
    if (allocateLabLocally)
        DYNALLOC2(int, labOutput, lab_sz, n, mWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, mWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, mWords, "malloc");

    stringtograph(tempg6, g, mWords); /// g6 string to densenauty
    GraphContainer resultContainer(n, mWords, g, rootedVertices.size()); /// container for graphs

    for (int i=0; i<rootedVertices.size(); ++i)
        resultContainer.SetRootedVertex(i, rootedVertices[i]);

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canong

    SetVertexColors(c, rootedVertices, n);
    SetColoredPartition(c, labOutput, ptn, n); /// set coloring (will be written over in call to densenauty)

    /// call densenauty
    densenauty(g, labOutput, ptn, orbits, &options, &stats, mWords, n, cg);

    if (stats.errstatus!=0)
        std::cerr << "AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty: call to NAUTY produced an error!\n";

    resultContainer.ColoredCanonicalRelabeling(labOutput, rootedVertices); /// relabel

    resultContainer.SetSymmFactor(AuxiliaryRoutinesForNauty::ExtractSymmFactor(stats)); /// get symmetry factor

#ifdef DEBUG
    GraphContainer refContainer(n, mWords, cg, rootedVertices.size());
    for (int i=0; i<rootedVertices.size(); ++i)
        refContainer.SetRootedVertex(i, i);
    if (refContainer==resultContainer)
        std::cout << "Success! resultContainer equals refContainer!\n";
    else
        std::cout << "ERROR: resultContainer does not equal refContainer!\n";
#endif

    free(c); /// free vertex colors

    delete[] tempg6;

    DYNFREE(g, g_sz); /// free graph
    DYNFREE(cg, cg_sz); /// free canonical graph
    if (allocateLabLocally) /// deallocate memory and set pointer back to NULL
    {
        DYNFREE(labOutput,lab_sz);
        labOutput = NULL;
    }
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    return resultContainer;
}

/// extract the symmetry factor from the statsblk structure
int AuxiliaryRoutinesForNauty::ExtractSymmFactor(statsblk status)
{
    /// write out to temporary file and read back in
    std::string tempFilename = "temp_symm"; /// HARD-CODED
    FILE *fpo = fopen(tempFilename.c_str(), "w"); /// open file for writing
    if (fpo==NULL)
        throw std::invalid_argument("ExtractSymmFactor: Error opening "+tempFilename);

    writegroupsize(fpo,status.grpsize1,status.grpsize2); /// write symmetry factor to file
    fprintf(fpo, "\n");

    fclose(fpo); /// close file

    FILE *fpi = fopen(tempFilename.c_str(), "r"); /// open file for reading
    if (fpi==NULL)
        throw std::invalid_argument("ExtractSymmFactor: Error opening "+tempFilename);

    int symmFactor = -1;
    fscanf(fpi, "%d\n", &symmFactor); /// read symmetry factor
    if (symmFactor==-1)
        throw std::invalid_argument("ExtractSymmFactor symmFactor has wrong value!\n");

    fclose(fpi); /// close file

    return symmFactor;
}
