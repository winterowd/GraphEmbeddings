#ifndef GRAPHGENERATORNAUTY_H
#define GRAPHGENERATORNAUTY_H

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <chrono>
#include <algorithm>

//#include "gtools.h"
#include "nauty.h"
#include "GraphContainer.h"

class GraphGeneratorParametersNauty /// graph generator parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// order of graph

    unsigned int L; /// number of bonds (FIXED)

    bool TwoRooted; /// do we generate two-rooted graphs from each connected graph?

    bool Verbose; /// verbosity

    bool ProcessCommandLine(int argc, char *argv[]); /// process command line arguments using BOOST

public:

    GraphGeneratorParametersNauty(int argc, char *argv[]); /// constructor

    /***** accessors *****/
    unsigned int GetN() const { return this->N; }

    unsigned int GetL() const { return this->L; }

    bool GenerateTwoRooted() const { return this->TwoRooted; }

    bool IsVerbose() const { return this->Verbose; }

};

class GraphGeneratorNauty
{
private:
    /***** private variables *****/
    GraphGeneratorParametersNauty Parameters;

    FILE* fp; /// file pointer for calling nauty routine readg

    int N; /// graph order

    int L; /// number of bonds (this is fixed)

    int MWords; /// nauty result from SETWORDSNEEDED(n)

    std::vector<std::vector<int>> RootedVertexNumbers; /// list of all possible rooted vertices where {v_1,v_2} \neq {v_2,v_1}

    /***** private routines *****/
    graph* GetNextGraph(graph *g); /// wrapper for readg

    void GetAllPossiblePairsForRootedVertices(); /// generate all possible pairs of rooted vertices

    void GenerateUniqueCombinationsWithNoDuplicates(std::vector<int>& tmp, const std::vector<int>& vertices, int k, bool verbose=false); /// generate combinations of integers contained in vertices with no repeats

    void GenerateTwoRootedFixedOrder(std::string inputFilename, bool verbose=false, bool outputSorted=false); /// generate rooted graphs from unrooted graphs of a given order

    void GenerateTwoRootedFixedOrderIterative(std::string inputFilename, bool verbose=false, bool outputSorted=false); /// more efficient algorithm which produces (hopefully) the same output as previous routine

    void ProduceNewLabelingGivenRootedVertices(const std::vector<int>& rooted, std::vector<int>& newLabeling, bool verbose=false); /// relabel vertices such that vertices of colors 0 and 1 are given labels 0 and 1, respectively

    void SetColoredPartition(int* c, int* lab, int* ptn); /// set up color parition of vertices

    void SetVertexColors(int *c, const std::vector<int>& rootedVertices, bool verbose=false); /// set the vertex colors for a given set of rooted vertices

    bool AreGraphParametersOK();

    /// debugging routine to compare to lists

public:
    GraphGeneratorNauty(int argc, char *argv[]); /// constructor

    void Generate(bool useIterativeRooted=true); /// interface for generating graphs

    void TestRelabeling(int n, std::string inputFilename); /// debugging routine
};

#endif // GRAPHGENERATORNAUTY_H
