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
#include "AuxiliaryRoutinesForNauty.h"

extern "C" {
int
GENG_MAIN(int argc, char *argv[]);
}

#include <boost/program_options.hpp>
namespace po = boost::program_options;

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

    bool AreGraphParametersOK();

    bool IsCanonical(graph *g); /// debugging routine to see if geng output is canonical

    /// debugging routine to compare to lists

public:

    int GetN() const { return this->N; }

    int GetL() const { return this->L; }

    bool GeneratedRooted() const { return this->Parameters.GenerateTwoRooted(); }

    GraphGeneratorNauty(int argc, char *argv[]); /// constructor

    void Generate(bool useIterativeRooted=true); /// interface for generating graphs

    std::vector<GraphContainer> GetColoredGraphsFromContainer(const GraphContainer& container); /// TODO: write this using GenerateTwoRootedFixedOrderIterative or GenerateTwoRootedFixedOrder

    /// get hardcoded filename for result of geng
    std::string GetOutputFilenameFixedOrder(bool rooted=false) const { if (rooted) return "graphs_g6_rooted_connected_N_"+std::to_string(this->N)+"_L_"+std::to_string(this->L)+".dat"; else return "graphs_g6_connected_N_"+std::to_string(this->N)+"_L_"+std::to_string(this->L)+".dat"; }
};

#endif // GRAPHGENERATORNAUTY_H
