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

#include "GraphContainer.h"
#include "AbstractLattice.h"
#include "SquareLattice.h"
#include "CubicLattice.h"
#include "TriangularLattice.h"

#include "gtools.h"
#include "nauty.h"

#include <iostream>

template<class Type> struct S; /// incomplete struct template for determiningn type of auto!

/// compute embedding number for simple connected undirected graphs on a given lattice

struct VertexEmbed {
    int Number; /// vertex number
    int Index; /// lattice index
};

inline std::ostream& operator<<(std::ostream& os, const VertexEmbed& v)
{
    os << v.Number << " at site " << v.Index;
    return os;
}

/// equality and inequality based on Index
inline bool operator==(const VertexEmbed& lhs, const VertexEmbed& rhs)
{
    return (lhs.Index == rhs.Index && lhs.Number == rhs.Number);
}

inline bool operator!=(const VertexEmbed& lhs, const VertexEmbed& rhs)
{
    return !(lhs==rhs);
}

inline bool operator==(const VertexEmbed& lhs, const int& rhs)
{
    return (lhs.Index == rhs);
}

inline bool operator!=(const VertexEmbed& lhs, const int& rhs)
{
    return !(lhs == rhs);
}

inline bool operator==(const int& lhs,  const VertexEmbed& rhs)
{
    return (lhs == rhs.Index);
}

inline bool operator!=(const int& lhs, const VertexEmbed& rhs)
{
    return !(lhs == rhs);
}

/// class to hold lists for embedding
/// when including next nearest neighbors and so on need to keep track of number of number of types of links used
/// getters and setters mostly. What other types of routines? any other data fields?
class VertexEmbedList
{
private:
    std::vector<VertexEmbed> List;

    std::vector<int> BondCounts; /// only need to keep track of total counts and not which edges are connected by NN or NNN

    int NbrChoicesForFirstBond; /// remember number of choices for first bond for unrooted graphs or 1 for rooted graphs

    bool TwoPointFunction; /// flag for two-point function (rooted graph)

     /// NOTE: when making canonical form, first vertex corresponds to color 1 and second vertex to color 2 (equality operator will distinguish them)
    std::vector<VertexEmbed> FixedVertices; /// for two-point functions (rooted graph)

    MaxInteractionLength CorrelatorDistance; /// keep track of which type of correlator

public:

    VertexEmbedList(MaxInteractionLength maxLength); /// constructor for unrooted graph

    VertexEmbedList(MaxInteractionLength maxLength, MaxInteractionLength correlatorDistance); /// constructor for correlator (rooted graph)

    VertexEmbedList(const VertexEmbedList& list) = default;

    void AddVertexEmbed(const VertexEmbed& v);

    void AddVertexEmbed(int number, int index);

    void AddFixedVerticesEmbed(const std::vector<VertexEmbed>& embed);

    void IncrementBondCount(int dIndex);

    int GetBondCount(int dIndex) const;

    VertexEmbed GetFixedVertex(int index) const;

    int GetSize() const { return this->List.size(); }

    int GetNbrBondTypes() const { return this->BondCounts.size(); }

    bool IsTwoPointFunction() const { return this->TwoPointFunction; }

    MaxInteractionLength GetCorrelatorDistance() const { return this->CorrelatorDistance; }

    int GetCorrelatorDistanceAsIndex() const { return static_cast<int>(this->CorrelatorDistance); }

    bool HasRepeatedVertices(); /// debugging routine

    bool HasRepeatedSites(); /// debugging routine

    VertexEmbed GetVertexEmbed(int index) const;

    void SetNbrChoicesForFirstBond(int nbr) { this->NbrChoicesForFirstBond = nbr; }

    int GetNbrChoicesForFirstBond() const { return this->NbrChoicesForFirstBond; }

    /// define iterator types
    using iterator = std::vector<VertexEmbed>::iterator;
    using const_iterator = std::vector<VertexEmbed>::const_iterator;

    /// begin and end functions
    iterator begin() { return this->List.begin(); }
    const_iterator begin() const { return this->List.begin(); }
    iterator end() { return this->List.end(); }
    const_iterator end() const { return this->List.end(); }

    friend std::ostream& operator<<(std::ostream& os, const VertexEmbedList& list); /// output for debugging reasons

};

bool operator==(const VertexEmbedList& lhs, const VertexEmbedList& rhs); /// comparison operator
bool operator!=(const VertexEmbedList& lhs, const VertexEmbedList& rhs); /// comparison operator

class GraphEmbedderParametersNauty /// graph embedder parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// order of graph

    std::string InputFilename; /// output of NAUTY geng for given order

    std::string OutputFilename; /// where to put embeddings and symmetry numbers?

    std::string InputFilenameFixedVertices; /// DEBUG: for now input these in a separate file of equal length

    std::string LatticeType; /// which lattice for embeddings

    MaxInteractionLength MaxEmbeddingLength; /// length of longest bond for embedding

    MaxInteractionLength CorrelatorLength; /// correlator length

    bool Correlator; /// flag for embedding correlator

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    GraphEmbedderParametersNauty(int argc, char *argv[]);

    unsigned int GetN() const { return this->N; }

    std::string GetInputFilename() const { return this->InputFilename; }

    std::string GetInputFilenameFixedVertices() const { return this->InputFilenameFixedVertices; }

    std::string GetOutputFilename() const { return this->OutputFilename; }

    std::string GetLatticeType() const { return this->LatticeType; }

    MaxInteractionLength GetMaxEmbeddingLength() const { return this->MaxEmbeddingLength; }

    MaxInteractionLength GetCorrelatorLength() const { return this->CorrelatorLength; }

    bool EmbedCorrelator() const { return this->Correlator; }

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

    enum { NbrLevelsNeighbor = 4 }; //// DEBUG: this is for array of function pointers! Change as capability for third nearest and further are added!

    int (GraphEmbedder::*GetNeighborFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    bool (GraphEmbedder::*AreNeighborsFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    std::array<int, NbrLevelsNeighbor> NbrNeighbors;

    std::vector<std::vector<int>> BondCombinations; /// for choosing combinations of links

    std::vector<std::vector<int>> FixedVertexNumbers; /// for choosing fixed vertices if embedding for two-point function (WRONG: need rooted graphs and not just trying all random pairs of a simple graph)
    /// NOTE: when making canonical form, first vertex corresponds to color 1 and second vertex to color 2

    /***** private methods *****/

    LatticeType ResolveLatticeType(const std::string& type); /// string to LatticeType

    graph* GetNextGraph(graph *g); /// wrapper for nauty routine readg

    int ComputeEmbeddingNumberNN(const GraphContainer& container); /// compute the embedding number of a graph where links are PURELY NN!

    void ComputeEmbeddingNumbers(const GraphContainer& container, FILE *fpo); /// loops over all valid combos of bond counts consistent with graph in container

    std::vector<int> ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo); /// compute the embedding number of a graph for a given combo of bond counts

    std::vector<int> GetAllowedBondDegreesOfList(const VertexEmbedList& list, const std::vector<int> &bondCombo); /// determine list of allowed neighbor degrees based on list and combination of bond counts

    bool IsProposedNNSiteFree(const std::vector<VertexEmbed>& list, int elem, int nn, int &newIndex);

    bool IsProposedNeighborSiteFree(const VertexEmbedList& list, int elem, int degree, int nn, int &newIndex);

    bool IsProposedSiteConsistentWithPreviousVerticesNN(const std::vector<VertexEmbed> &list, const GraphContainer &container, int elem, int newIndex, int newVertexNumber);

    bool IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(const VertexEmbedList& list, const GraphContainer &container, const std::vector<int> &bondCombo, int elem, int newIndex, int newVertexNumber, std::vector<int> &bondCountsToBeAdded);

    std::vector<VertexEmbed> SelectFirstLinkToEmbedNN(const GraphContainer& container);

    std::vector<VertexEmbedList> CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::vector<VertexEmbedList> CreateInitialVertexEmbedListsRootedFixed(const GraphContainer& container, const std::vector<int> &bondCombo, const std::vector<int> &fixedVertices);

    std::vector<VertexEmbedList> CreateInitialVertexEmbedListsRooted(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::vector<VertexEmbedList> CreateInitialVertexEmbedListsNonRooted(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::unordered_set<int> GetRemainingVertices(const std::vector<VertexEmbed>& listUsedVertices);

    std::unordered_set<int> GetRemainingVertices(const VertexEmbedList& listUsedVertices);

    void GetCombinationsOfBonds(int nbrBonds); /// compute numbers of each type of bonds given the total number of bonds for a graph which is to be embedded

    void GenerateCombinations(const std::vector<int>& arr, std::vector<int>& data, int index); /// called by GetCombinationsOfBonds

    void GetAllPossiblePairsForFixedVertices(); /// all choices of fixed vertices for two point function (ASSUMES: NN)

    void GenerateUniqueCombinationsWithNoDuplicates(std::vector<int>& tmp, const std::vector<int>& vertices, int left, int k);

    bool IsDuplicate(const std::vector<std::vector<VertexEmbed>>& lists, const std::vector<VertexEmbed>& toAdd);

    bool AreBondCountsEqual(const VertexEmbedList& lhs, const VertexEmbedList& rhs);

    bool IsDuplicateOld(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose=false); /// will use in embedding routine

    bool IsDuplicate(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose=false);

    void TestEraseWrongSizes(std::vector<VertexEmbedList>& lists, int vertexCount); /// debugging routine

    void CallDenseNauty(graph *g, int *lab, int *ptn, int *orbits, statsblk &stats);

    std::vector<int> GetFixedVertexPair(int index) { if (index>=0 && index<=this->FixedVertexNumbers.size()) return this->FixedVertexNumbers[index]; return std::vector<int>{-1}; }

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

public:

    GraphEmbedder(int argc, char *argv[]); /// new style constructor

    ~GraphEmbedder(); /// destructor

    void Embed(bool debugJonas=false); /// calculates embedding numbers for all graphs of a given order and outputs results

    void EmbedSpecificGraphBondCombo(int graphNbr, const std::vector<int>& bondCounts); /// debugging routine: graphNbr is the line in the file and bondCounts needs to be provided by the user

    void TestInitialRootedGraphList(std::string filename); /// test routine that reads in a graph in g6 format in filename and fixes v1 and v2 as rooted

};

#endif // GRAPHEMBEDDER_H
