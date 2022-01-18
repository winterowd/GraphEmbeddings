#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

extern "C" {
#include "nauty.h"
}
#include "gtools.h"
#include "AuxiliaryRoutinesForNauty.h"
#include "MiscellaneousRoutines.h"

/// read in Jonas' output file (all graphs with six bonds) and reformat it
/// @param inputFilename: name of file
/// @param nMax: maximum number of vertices
/// INPUT: first line (ignore)
/// g6_string rooted_vertices (0,6) (1,5) (2,4) (3,3) (4,2) (5,1) (6,0)
/// OUTPUT: canonical_g6_string (0,6) (1,5) (2,4) (3,3) (4,2) (5,1) (6,0)
void ReadAndReformatJonas(std::string inputFilename, int nMax=10)
{
    std::ifstream jonasFile(inputFilename);

    std::string tempLine;

    std::getline(jonasFile, tempLine); /// ignore first line

    int m = SETWORDSNEEDED(nMax);

    nauty_check(WORDSIZE, nMax, m, NAUTYVERSIONID);

    /// allocate and declare NAUTY stuff
    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graphs
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    DYNALLOC2(graph, g, g_sz, nMax, m, "malloc");
    DYNALLOC2(graph, cg, cg_sz, nMax, m, "malloc");
    DYNALLOC2(int, lab, lab_sz, nMax, m, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, nMax, m, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, nMax, m, "malloc");

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canong

    int count=2;
    while (std::getline(jonasFile, tempLine))
    {
        std::string g6String;
        std::vector<int> tempRootedList;
        std::string token;

        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string

        /// convert std::string to c-style string
        char *tempg6 = new char[g6String.length()+1];
        std::strcpy(tempg6, g6String.c_str());

        int n = graphsize(tempg6); /// get current size
        stringtograph(tempg6, g, m); /// g6 string to densenauty

        delete[] tempg6; /// clean up memory

        /// read in rooted vertices (start at 1 in Jonas' convention instead of 0 in our convention)
        ss >> token;
        token.erase(std::remove(token.begin(), token.end(), '{'), token.end());
        token.erase(std::remove(token.begin(), token.end(), ','), token.end());
        tempRootedList.push_back(std::stoi(token)-1);

        ss >> token;
        token.erase(std::remove(token.begin(), token.end(), '}'), token.end());
        tempRootedList.push_back(std::stoi(token)-1);

        int *c = (int*)malloc(n * sizeof(int)); /// alloc C-style arrray for colors of vertices

        AuxiliaryRoutinesForNauty::SetVertexColors(c, tempRootedList, n); /// set vertex colors based on rooted vertices
        AuxiliaryRoutinesForNauty::SetColoredPartition(c, lab, ptn, n); /// set up lab and ptn based on colors

        densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg); /// call nauty

        std::string g6canonicalString(ntog6(cg,m,n)); /// convert canonical graph to std::string
        g6canonicalString.erase(std::remove(g6canonicalString.begin(), g6canonicalString.end(), '\n'), g6canonicalString.end()); // get rid of newline character

        free (c); /// clean up

        std::vector<int> embeddings; /// read in embeddings (do nothing with them)
        for (int x; ss >> x;)
            embeddings.push_back(x);

        if (embeddings.size() != 7) /// partitions of 6 bonds into NN and NNN
            throw std::invalid_argument("SHOULD HAVE 7 PARTITIONS!\n");

        std::cout << g6canonicalString;
        for (auto x: embeddings)
            std::cout << " " << x;
        std::cout << "\n";

    }

    DYNFREE(g, g_sz);
    DYNFREE(cg, cg_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    jonasFile.close();

}

int main(int argc, char *argv[])
{
    ReadAndReformatJonas("jonas_embeddings.dat");
    return 0;
}
