#include <iostream>
#include <vector>

#include "AuxiliaryRoutinesForNauty.h"
extern "C" {
#include "nauty.h"
}
#include "gtools.h"

/// test routine for dense nauty canonicalization on rooted graphs
/// takes the fifth-order graph "DrO" and chooses vertices (1,4) (4,1) (2,0) and (2,3) to be rooted
/// after canonicalization of each rooted graph, the first two choices for the rooted vertices should produce different canonical graphs and the second two choices should produce the same canonical graph
/// first rooted vertex should be relabeled as vertex 0 and second rooted vertex should be relabeled as vertex 1
void TestDenseNautyColoredCanon(bool verbose=false)
{
    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(graph, cg1, cg1_sz);
    DYNALLSTAT(graph, cg2, cg2_sz);
    DYNALLSTAT(graph, cg3, cg3_sz);
    DYNALLSTAT(graph, cg4, cg4_sz);
    DYNALLSTAT(int, lab1, lab1_sz);
    DYNALLSTAT(int, lab2, lab2_sz);
    DYNALLSTAT(int, lab3, lab3_sz);
    DYNALLSTAT(int, lab4, lab4_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    statsblk stats;

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canonical

    int n = 5;
    int m = SETWORDSNEEDED(n); /// array of m setwords sufficients to hold n bits

    nauty_check(WORDSIZE, n, m, NAUTYVERSIONID); /// check if everything is ok

    DYNALLOC2(graph, g, g_sz, n, m, "malloc");
    DYNALLOC2(graph, cg1, cg1_sz, n, m, "malloc");
    DYNALLOC2(graph, cg2, cg2_sz, n, m, "malloc");
    DYNALLOC2(graph, cg3, cg3_sz, n, m, "malloc");
    DYNALLOC2(graph, cg4, cg4_sz, n, m, "malloc");
    DYNALLOC2(int, lab1, lab1_sz, n, m, "malloc");
    DYNALLOC2(int, lab2, lab2_sz, n, m, "malloc");
    DYNALLOC2(int, lab3, lab3_sz, n, m, "malloc");
    DYNALLOC2(int, lab4, lab4_sz, n, m, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, m, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, m, "malloc");

    std::vector<int> tempRootedList{-1,-1}; /// only two rooted vertices

    int *c = (int*)malloc(n * sizeof(int)); /// alloc C-style arrray for colors of vertices

    /// set vertex colors based on rooted vertices and then set up lab and ptn

    /// first two should give different canonical graphs
    tempRootedList[0] = 1; tempRootedList[1] = 4;
    SetVertexColors(c, tempRootedList, n);
    SetColoredPartition(c, lab1, ptn, n);

    tempRootedList[0] = 4; tempRootedList[1] = 1;
    SetVertexColors(c, tempRootedList, n);
    SetColoredPartition(c, lab2, ptn, n);

    /// second two should give identical canonical graphs
    tempRootedList[0] = 2; tempRootedList[1] = 0;
    SetVertexColors(c, tempRootedList, n);
    SetColoredPartition(c, lab3, ptn, n);

    tempRootedList[0] = 2; tempRootedList[1] = 3;
    SetVertexColors(c, tempRootedList, n);
    SetColoredPartition(c, lab4, ptn, n);

    free(c); /// free c-style array

    EMPTYGRAPH(g, n, m); /// clear graph

    /// set up graph
    ADDONEEDGE(g, 0, 1, m);
    ADDONEEDGE(g, 0, 2, m);
    ADDONEEDGE(g, 2, 3, m);
    ADDONEEDGE(g, 3, 1, m);
    ADDONEEDGE(g, 1, 4, m);

    set *gj; /// for output

    if (verbose)
    {
        for (int j = 1; j < n; ++j)
        {
            gj = GRAPHROW(g,j,m);
            for (int i = 0; i < j; ++i)
            {
                if (ISELEMENT(gj,i))
                    std::cout << "G: Vertex " << j << " connected to vertex " << i << "\n";
            }
        }
    }

    densenauty(g, lab1, ptn, orbits, &options, &stats, m, n, cg1);
    densenauty(g, lab2, ptn, orbits, &options, &stats, m, n, cg2);
    densenauty(g, lab3, ptn, orbits, &options, &stats, m, n, cg3);
    densenauty(g, lab4, ptn, orbits, &options, &stats, m, n, cg4);

    for (int i=0; i<n; ++i)
        std::cout << "lab1: vertex " << lab1[i] << " maps to vertex " << i << "\n";

    for (int i=0; i<n; ++i)
        std::cout << "lab2: vertex " << lab2[i] << " maps to vertex " << i << "\n";

    for (int i=0; i<n; ++i)
        std::cout << "lab3: vertex " << lab3[i] << " maps to vertex " << i << "\n";

    for (int i=0; i<n; ++i)
        std::cout << "lab4: vertex " << lab4[i] << " maps to vertex " << i << "\n";

    if (verbose)
    {
        for (int j = 1; j < n; ++j)
        {
            gj = GRAPHROW(cg1,j,m);
            for (int i = 0; i < j; ++i)
            {
                if (ISELEMENT(gj,i))
                    std::cout << "CG1: Vertex " << j << " connected to vertex " << i << "\n";
            }
        }
        for (int j = 1; j < n; ++j)
        {
            gj = GRAPHROW(cg2,j,m);
            for (int i = 0; i < j; ++i)
            {
                if (ISELEMENT(gj,i))
                    std::cout << "CG2: Vertex " << j << " connected to vertex " << i << "\n";
            }
        }
        for (int j = 1; j < n; ++j)
        {
            gj = GRAPHROW(cg3,j,m);
            for (int i = 0; i < j; ++i)
            {
                if (ISELEMENT(gj,i))
                    std::cout << "CG3: Vertex " << j << " connected to vertex " << i << "\n";
            }
        }
        for (int j = 1; j < n; ++j)
        {
            gj = GRAPHROW(cg4,j,m);
            for (int i = 0; i < j; ++i)
            {
                if (ISELEMENT(gj,i))
                    std::cout << "CG4: Vertex " << j << " connected to vertex " << i << "\n";
            }
        }
    }

    DYNFREE(lab1,lab1_sz);
    DYNFREE(lab2,lab2_sz);
    DYNFREE(lab3,lab3_sz);
    DYNFREE(lab4,lab4_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);
    DYNFREE(g,g_sz);
    DYNFREE(cg1,cg1_sz);
    DYNFREE(cg2,cg2_sz);
    DYNFREE(cg3,cg3_sz);
    DYNFREE(cg4,cg4_sz);

}

int main()
{
    std::cout << "RUNNING: dense NAUTY test on 5th order graph DrO!\n";
    TestDenseNautyColoredCanon(true);
    return 0;
}
