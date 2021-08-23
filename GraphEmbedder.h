#ifndef GRAPHEMBEDDER_H
#define GRAPHEMBEDDER_H

#include <string>
#include <vector>
#include <unordered_set>
#include <array>
#include <numeric>
#include <chrono>

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

    int NbrChoicesForFirstBond; /// remember type of first bond

public:

    VertexEmbedList(MaxInteractionLength maxLength);

    VertexEmbedList(const VertexEmbedList& list) = default;

    void AddVertexEmbed(const VertexEmbed& v);

    void AddVertexEmbed(int number, int index);

    void IncrementBondCount(int dIndex);

    int GetBondCount(int dIndex) const;

    int GetSize() const { return this->List.size(); }

    int GetNbrBondTypes() const { return this->BondCounts.size(); }

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

class GraphEmbedderParametersNauty /// graph embedder parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// order of graph

    std::string InputFilename; /// output of NAUTY geng for given order

    std::string OutputFilename; /// where to put embeddings and symmetry numbers?

    std::string LatticeType; /// which lattice for embeddings

    MaxInteractionLength Length; /// length of longest bond for embedding

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    GraphEmbedderParametersNauty(int argc, char *argv[]);

    unsigned int GetN() const { return this->N; }

    std::string GetInputFilename() const { return this->InputFilename; }

    std::string GetOutputFilename() const { return this->OutputFilename; }

    std::string GetLatticeType() const { return this->LatticeType; }

    MaxInteractionLength GetMaxInteractionLength() const { return this->Length; }

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

    enum { NbrLevelsNeighbor = 2 }; //// DEBUG: this is for array of function pointers! Change is third nearest and further are added!

    int (GraphEmbedder::*GetNeighborFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    bool (GraphEmbedder::*AreNeighborsFunctionPointerArray[NbrLevelsNeighbor])(int, int);

    std::array<int, NbrLevelsNeighbor> NbrNeighbors;

    std::vector<std::vector<int>> BondCombinations; /// for choosing combinations of links

    /***** private methods *****/

    LatticeType ResolveLatticeType(const std::string& type); /// string to LatticeType

    graph* GetNextGraph(graph *g); /// wrapper for nauty routine readg

    int ComputeEmbeddingNumberNN(const GraphContainer& container); /// compute the embedding number of a graph where links are PURELY NN!

    void ComputeEmbeddingNumbers(const GraphContainer& container, FILE *fpo); /// loops over all valid combos of bond counts consistent with graph in container

    int ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo); /// compute the embedding number of a graph for a given combo of bond counts

    std::vector<int> GetAllowedBondDegreesOfList(const VertexEmbedList& list, const std::vector<int> &bondCombo); /// determine list of allowed neighbor degrees based on list and combination of bond counts

    bool IsProposedNNSiteFree(const std::vector<VertexEmbed>& list, int elem, int nn, int &newIndex);

    bool IsProposedNeighborSiteFree(const VertexEmbedList& list, int elem, int degree, int nn, int &newIndex);

    bool IsProposedSiteConsistentWithPreviousVerticesNN(const std::vector<VertexEmbed> &list, const GraphContainer &container, int elem, int newIndex, int newVertexNumber);

    bool IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(const VertexEmbedList& list, const GraphContainer &container, const std::vector<int> &bondCombo, int elem, int newIndex, int newVertexNumber, std::vector<int> &bondCountsToBeAdded);

    std::vector<VertexEmbed> SelectFirstLinkToEmbedNN(const GraphContainer& container);

    std::vector<VertexEmbedList> CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo);

    std::unordered_set<int> GetRemainingVertices(const std::vector<VertexEmbed>& listUsedVertices);

    std::unordered_set<int> GetRemainingVertices(const VertexEmbedList& listUsedVertices);

    void GetCombinationsOfBonds(int nbrBonds); /// compute numbers of each type of bonds given the total number of bonds for a graph which is to be embedded

    void GenerateCombinations(const std::vector<int>& arr, std::vector<int>& data, int index); /// called by GetCombinationsOfBonds

    bool IsDuplicate(const std::vector<std::vector<VertexEmbed>>& lists, const std::vector<VertexEmbed>& toAdd);

    bool AreBondCountsEqual(const VertexEmbedList& lhs, const VertexEmbedList& rhs);

    bool IsDuplicate(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose=false); /// will use in embedding routine

    void TestEraseWrongSizes(std::vector<VertexEmbedList>& lists, int vertexCount); /// debugging routine

    void CallDenseNauty(graph *g, int *lab, int *ptn, int *orbits, statsblk &stats);

    /// wrapper for lattice accessor
    int GetNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetNearestNeighbor(siteIndex, nnIndex); }

    /// wrapper for lattice accessor
    int GetNextNearestNeighbor(int siteIndex, int nnIndex) { return this->Lattice->GetNextNearestNeighbor(siteIndex, nnIndex); }

     /// wrapper for lattice accessor
    bool AreNN(int index1, int index2) { return this->Lattice->AreNN(index1, index2); }

    /// wrapper for lattice accessor
    bool AreNNN(int index1, int index2) { return this->Lattice->AreNNN(index1, index2); }

    int GetNeighbor(int degree, int siteIndex, int neighborIndex);

    int GetNbrNeighbors(int degree);

    bool AreNeighbors(int degree, int index1, int index2);

    void UpdateBondCounts(VertexEmbedList &List, const std::vector<int>& countsToBeAdded);

public:

    GraphEmbedder(int argc, char *argv[]); /// new style constructor

    ~GraphEmbedder(); /// destructor

    void Embed(); /// calculates embedding numbers for all graphs of a given order and outputs results 

    void EmbedSpecificGraphBondCombo(int graphNbr, const std::vector<int>& bondCounts); /// debugging routine: graphNbr is the line in the file and bondCounts needs to be provided by the user

};

#endif // GRAPHEMBEDDER_H
