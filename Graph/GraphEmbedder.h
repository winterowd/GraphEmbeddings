#ifndef GRAPHEMBEDDER_H
#define GRAPHEMBEDDER_H

#include <string>
#include <vector>
#include <unordered_set>
#include <array>
#include <numeric>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tuple>

/// user-defined headers
#include "GraphContainer.h"
#include "AbstractLattice.h"
#include "SquareLattice.h"
#include "CubicLattice.h"
#include "TriangularLattice.h"
#include "VertexEmbedList.h"
#include "CubicLatticeCanonicalizor.h"

/// nauty headers
#include "gtools.h"
#include "nauty.h"

#define MAX_G6_LENGTH 50

template<class Type> struct S; /// incomplete struct template for determiningn type of auto!

/// compute embedding number for simple connected undirected graphs on a given lattice

class GraphEmbedderParametersNauty /// graph embedder parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// order of graph

    std::string InputFilename; /// output of NAUTY geng for given order

    std::string OutputFilename; /// where to put embeddings and symmetry numbers?

    std::string LatticeType; /// which lattice for embeddings

    std::string G6; /// g6 string

    MaxInteractionLength MaxEmbeddingLength; /// length of longest bond for embedding

    MaxInteractionLength CorrelatorLength; /// correlator length

    bool Correlator; /// flag for embedding correlator

    bool JonasFormat; /// output embeddings in Jonas' format

    void RequiredOptionWhenOtherOptionMissing(const po::variables_map& vm, const char* required, const char* dependancy);

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    GraphEmbedderParametersNauty(int argc, char *argv[]);

    unsigned int GetN() const { return this->N; }

    std::string GetInputFilename() const { return this->InputFilename; }

    std::string GetOutputFilename() const { return this->OutputFilename; }

    std::string GetLatticeType() const { return this->LatticeType; }

    std::string GetG6() const { return this->G6; }

    MaxInteractionLength GetMaxEmbeddingLength() const { return this->MaxEmbeddingLength; }

    MaxInteractionLength GetCorrelatorLength() const { return this->CorrelatorLength; }

    bool EmbedCorrelator() const { return this->Correlator; }

    bool UseJonasFormat() const { return this->JonasFormat; }

};


class GraphEmbedder
{
private:

    enum LatticeType {
        Square,
        Triangular,
        Cubic,
        Invalid
    };

    /***** private variables *****/
    GraphEmbedderParametersNauty Parameters;

    int N; /// number of vertices

    int MWords; /// nauty result from SETWORDSNEEDED(n)

    int MaxDegreeNeighbor; /// maximum of degree of neighbors for embedding

    FILE* fp; /// file pointer for calling nauty routine readg

    AbstractLattice *Lattice; /// lattice object

    std::unordered_set<int> VertexSet; /// precompute for embedding number

    bool DebugJonas; /// debug flag

    std::string G6String; /// string if we want to access mode where a single graph is embedded by inputting the g6 string

    enum { NbrLevelsNeighbor = 4 }; //// DEBUG: this is for array of function pointers! Change as capability for third nearest and further are added!

    int (GraphEmbedder::*GetNeighborFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    bool (GraphEmbedder::*AreNeighborsFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    std::array<int, NbrLevelsNeighbor> NbrNeighbors;

    std::vector<std::vector<int>> BondCombinations; /// for choosing combinations of links

    std::set<VertexEmbedList> EmbedLists; /// this should get set each time we call ComputeEmbeddingNumberCombo and contain the lists

    /// NOTE: when making canonical form, colors 0 and 1 correspond to vertices 1 and 2, respectively
    ///
    /***** private methods *****/

    LatticeType ResolveLatticeType(const std::string& type); /// string to LatticeType

    graph* GetNextGraph(graph *g); /// wrapper for nauty routine readg

    int ComputeEmbeddingNumberNN(const GraphContainer& container); /// compute the embedding number of a graph where links are PURELY NN!

    void ComputeEmbeddingNumbers(const GraphContainer& container, graph *g, FILE *fpo, int symmFactor=-1); /// loops over all valid combos of bond counts consistent with graph in container

    int ComputeEmbeddingNumberComboOldWorking(const GraphContainer& container, const std::vector<int> &bondCombo); /// compute the embedding number of a graph for a given combo of bond counts

    std::pair<int,VertexEmbedList> ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo); /// compute the embedding number of a graph for a given combo of bond counts

    std::pair<std::vector<VertexEmbedList>, std::vector<int>> ComputeCanonicalGraphsAndEmbeddingNumbers(GraphContainer container, const std::vector<int> &bondCombo);

    std::vector<int> GetAllowedBondDegreesOfList(const VertexEmbedList& list, const std::vector<int> &bondCombo); /// determine list of allowed neighbor degrees based on list and combination of bond counts

    bool IsProposedNNSiteFree(const std::vector<VertexEmbed>& list, int elem, int nn, int &newIndex);

    bool IsProposedNeighborSiteFree(const VertexEmbedList& list, VertexEmbed elem, int degree, int nn, int &newIndex);

    bool IsProposedSiteConsistentWithPreviousVerticesNN(const std::vector<VertexEmbed> &list, const GraphContainer &container, int elem, int newIndex, int newVertexNumber);

    bool IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(const VertexEmbedList& list, const GraphContainer &container, const std::vector<int> &bondCombo, VertexEmbed elem, int newIndex, int newVertexNumber, std::vector<int> &bondCountsToBeAdded);

    std::vector<VertexEmbed> SelectFirstLinkToEmbedNN(const GraphContainer& container);

    std::set<VertexEmbedList> CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::set<VertexEmbedList> CreateInitialVertexEmbedListsRootedFixed(const GraphContainer& container, const std::vector<int> &bondCombo, const std::vector<int> &rootedVertices);

    std::set<VertexEmbedList> CreateInitialVertexEmbedListsRooted(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::set<VertexEmbedList> CreateInitialVertexEmbedListsNonRooted(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::unordered_set<int> GetRemainingVertices(const std::vector<VertexEmbed>& listUsedVertices);

    std::unordered_set<int> GetRemainingVertices(const VertexEmbedList& listUsedVertices);

    void GetCombinationsOfBondsFixedManhattanDistance(int nbrBonds, int manhattanDistance); /// compute numbers of each type of bonds with a constraint for a given "Manhattan Distance"

    void GetCombinationsOfBondsFixedNumberOfBonds(int nbrBonds); /// compute numbers of each type of bonds given the total number of bonds for a graph which is to be embedded

    void GenerateCombinations(const std::vector<int>& arr, std::vector<int>& data, int index, std::function<bool(const std::vector<int>&)> isValid);

    bool IsDuplicate(const std::vector<std::vector<VertexEmbed>>& lists, const std::vector<VertexEmbed>& toAdd);

    bool AreBondCountsEqual(const VertexEmbedList& lhs, const VertexEmbedList& rhs);

    bool IsDuplicateOld(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose=false); /// will use in embedding routine

    bool IsDuplicate(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd);

    bool IsDuplicate(const std::set<VertexEmbedList>& lists, const VertexEmbedList& toAdd);

    void TestEraseWrongSizes(std::vector<VertexEmbedList>& lists, int vertexCount); /// debugging routine

    std::pair<int, VertexEmbed> DetermineNextVertexToEmbed(const GraphContainer& container, const VertexEmbedList& embedded, const std::unordered_set<int>& remainingVertices);

    void CallDenseNauty(graph *g, int *lab, int *ptn, int *orbits, statsblk &stats);

    /// wrapper for lattice accessor
    int GetNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetNearestNeighbor(siteIndex, nnIndex); }

    /// wrapper for lattice accessor
    int GetNextNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetNextNearestNeighbor(siteIndex, nnIndex); }

    /// wrapper for lattice accessor
    int GetThirdNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetThirdNearestNeighbor(siteIndex, nnIndex); }

    /// wrapper for lattice accessor
    int GetFourthNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetFourthNearestNeighbor(siteIndex, nnIndex); }

    /// wrapper for lattice accessor
    bool AreNN(int index1, int index2) { return this->Lattice->AreNN(index1, index2); }

    /// wrapper for lattice accessor
    bool AreNNN(int index1, int index2) { return this->Lattice->AreNNN(index1, index2); }

    /// wrapper for lattice accessor
    bool AreThirdNN(int index1, int index2) { return this->Lattice->AreThirdNN(index1, index2); }

    /// wrapper for lattice accessor
    bool AreFourthNN(int index1, int index2) { return this->Lattice->AreFourthNN(index1, index2); }

    int GetNeighbor(int degree, int siteIndex, int neighborIndex);

    int GetNbrNeighbors(int degree);

    bool AreNeighbors(int degree, int index1, int index2);

    void UpdateBondCounts(VertexEmbedList &List, const std::vector<int>& countsToBeAdded);

    void PrintVertexEmbedList(const VertexEmbedList& list);

    void EmbedFromFile();

    void EmbedSingleG6();

    int GetGraphSizeFromString(const std::string& g6String);

    void DenseNautyFromString(const std::string& g6String, graph *g);

    int GetSymmFactor(graph *g);

public:

    GraphEmbedder(int argc, char *argv[]); /// new style constructor

    ~GraphEmbedder(); /// destructor

    void Embed(); /// calculates embedding numbers for all graphs of a given order or a single and outputs results

    /// TODO: routine which returns all partitions for a given container

    /// TODO: routine which takes a single partition and container returns all VertexEmbedLists

    /// routine which outputs vector of canonical graphs and counts for embeddings (all NN links)
    std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> GetCanonicalGraphsAndCounts(); /// uses the graph provided as a g6 string by the user (command line arguments sent to constructor)

    std::pair<GraphContainer, VertexEmbedList> ContainerAndSampleCubicEmbeddingFromG6(); /// return a pair consisting of a container and a VertexEmbedList for a given g6 string input by the user (for debugging purposes)

};

#endif // GRAPHEMBEDDER_H
