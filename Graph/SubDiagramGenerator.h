#ifndef SUBDIAGRAMGENERATOR_H
#define SUBDIAGRAMGENERATOR_H

#include <algorithm>

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "CubicLattice.h"
#include "CubicLatticeCanonicalizor.h"
#include "AuxiliaryRoutinesForNauty.h"

/// helper class for D_n: set of sets of n elements where each element is a connected subgraph of G and all are pairwise disjoint
template<class T>
class SetOfSets
{
private:
    int N; /// required size of each subset

    std::set<std::set<T>> Sets; /// set of sets
public:
    SetOfSets(int n) : N(n) {}; /// constructor

    void AddSet(std::set<T> set) /// add a set of size N
    {
        if (set.size()==this->N)
        {
            auto result = Sets.insert(set);
#ifdef DEBUG
            if (!result.second)
                std::cout << "AddSet: failed to insert!\n";
#endif
        }
        else
            std::cout << "ERROR: Set not of correct size()! " << set.size() << " " << this->N << "\n";
    }

    int GetN() const { return this->N; } /// get N

    int GetNbrSets() const { return this->Sets.size(); } /// get the number of sets

    /// define iterator types
    using iterator = typename std::set<std::set<T>>::iterator;
    using const_iterator = typename std::set<std::set<T>>::const_iterator;

    /// begin and end functions
    iterator begin() { return this->Sets.begin(); }
    const_iterator begin() const { return this->Sets.begin(); }
    iterator end() { return this->Sets.end(); }
    const_iterator end() const { return this->Sets.end(); }

};

template<class T>
inline bool operator<(const SetOfSets<T>& lhs,  const SetOfSets<T>& rhs)
{
    return (lhs.GetN()<rhs.GetN());
}

/// data structure for an embedded undirected edge
struct UndirectedEmbeddedEdge {
    int index1; /// first lattice index
    int index2; /// second lattice index
};

/// output edge
inline std::ostream& operator<<(std::ostream& os, const UndirectedEmbeddedEdge& e)
{
    os << "(" << e.index1 << "," << e.index2 << ")";
    return os;
}

/// equality and inequality operators for embedded edge
inline bool operator==(const UndirectedEmbeddedEdge& lhs, const UndirectedEmbeddedEdge& rhs)
{
    return ((lhs.index1 == rhs.index1 && lhs.index2 == rhs.index2) || (lhs.index1 == rhs.index2 && lhs.index2 == rhs.index1));
}

inline bool operator!=(const UndirectedEmbeddedEdge& lhs, const UndirectedEmbeddedEdge& rhs)
{
    return !(lhs==rhs);
}

/// equality operator for two sets of embedded edges
inline bool operator==(const std::vector<UndirectedEmbeddedEdge>& lhs, const std::vector<UndirectedEmbeddedEdge>& rhs)
{
    if (lhs.size()!=rhs.size())
        return false;
    for (auto it=lhs.begin(); it!=lhs.end(); ++it)
        if (std::find(rhs.begin(), rhs.end(), *it) == rhs.end())
            return false;
    return true;
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<UndirectedEmbeddedEdge>& edges)
{
    for (auto it=edges.begin(); it!=edges.end(); ++it)
        std::cout << *it << " ";
    return os;
}

/// equality for a vectors of edges (representing some embedded graphs)
inline bool operator!=(const std::vector<UndirectedEmbeddedEdge>& lhs, const std::vector<UndirectedEmbeddedEdge>& rhs)
{
    return !(lhs==rhs);
}

/// helper class in order to distinguish canonical subgraphs
class CanonicalSubDiagram
{
private:
    VertexEmbedList CanonicalList; /// canonical vertex embed list (output from CubicLatticeCanonicalizor)

    GraphContainer CanonicalContainer; /// canonical container (output from NAUTY)
public:
    CanonicalSubDiagram(const VertexEmbedList& list, const GraphContainer& container) : CanonicalList(list), CanonicalContainer(container) {}

    int GetN() const { return this->CanonicalContainer.GetN(); }

    int GetL() const { return this->CanonicalContainer.GetL(); }

    GraphContainer GetCanonicalContainer() const { return this->CanonicalContainer; }

    VertexEmbedList GetCanonicalList() const { return this->CanonicalList; }

    friend bool operator==(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs);
    friend std::ostream& operator<< (std::ostream& stream, const CanonicalSubDiagram& can);
};

/// equality operator (compare embed list and container!)
inline bool operator==(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs)
{
    return (lhs.CanonicalList==rhs.CanonicalList && lhs.CanonicalContainer==rhs.CanonicalContainer);
}

inline bool operator!=(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs)
{
   return !(lhs==rhs);
}

inline std::ostream& operator<< (std::ostream& stream, const CanonicalSubDiagram& can)
{
    stream << "CANONICAL_GRAPH:\n";
    stream << can.CanonicalContainer;
    stream << can.CanonicalList;
    return stream;
}

/// generate all subgraphs of a given graph (embedded on a lattice i.e. cubic)
class SubDiagramGenerator
{
private:
    /***** private variables *****/
    GraphContainer OriginalContainer; /// pointer to original graph container TODO: change this so that it has a copy

    VertexEmbedList OriginalList; /// pointer to original vertex embed list TODO: changes this so that it has a copy

    CubicLattice *MyCubicLattice; /// pointer to cubic lattice object

    std::vector<std::vector<std::pair<int, GraphContainer>>> SortedSubDiagramsWithMap; /// subdiagrams sorted by number of bonds together with (UNSORTED) index matching the vertex map

    /// ith member of VerticesMap contains an array of length N_{v,sg_i}
    /// where VerticesMap[i][k] contains the the ORIGINAL vertex label of the vertex RELABELED k+1 in the ith UNSORTED subgraph (k=0,1,...,N_{v,sg_i}-1)
    std::vector<std::vector<int>> VerticesMap; /// map of vertices of subgraphs (order differs with SortedSubDiagrams i.e. UNSORTED)

    std::vector<VertexEmbedList> EmbedLists; ///  embed list for the subgraphs (same order as VerticesMap i.e. UNSORTED)

    std::vector<VertexEmbedList> EmbedListsCanonicalLabels; /// embed list for the subgraphs with canonical vertex labels (SORTED)

    std::vector<std::vector<UndirectedEmbeddedEdge>> EmbeddedEdgeLists; /// edge lits for the subgraphs (UNSORTED)

    std::vector<int> SubgraphToCanonicalMap; /// element i contains index of element of CanonicalSubDiagramList to which subdiagram i corresponds to (SORTED)

    //// ith member of CanonicalToOriginalMap contains an array of of length N_{v,sg_i}
    /// where CanonicalToOriginalMap[i][k] contains the the ORIGINAL vertex label of the CANONICALLY labeled vertex k+1 in the ith SORTED subgraph (k=0,1,...,N_{v,sg_i}-1)
    std::vector<std::vector<int>> CanonicalToOriginalMap; /// mapping between canonical and original labels

    /// objects which contain container and lists which are canonical with respect to NAUTY and the octohedral group! (cubic symmetry)
    std::vector<CanonicalSubDiagram> CanonicalSubDiagramList;

    std::vector<SetOfSets<int>> DisjointSets; /// D_n: set of subsets represented by the index corresponding to ordering of subgraphs in SortedSubDiagramsWithMap
    int NbrSubDiagrams; /// total number of connected subdiagrams

    /**** private methods ****/

    //// canonicalize with respect to NAUTY and cubic symmetries
    std::pair<CanonicalSubDiagram, VertexEmbedList> ComputeCanonicalSubgraphAndListUnrooted(int sortedIndex);

    std::pair<CanonicalSubDiagram, VertexEmbedList> ComputeCanonicalSubgraphAndListRooted(int sortedIndex);

    std::pair<CanonicalSubDiagram, VertexEmbedList> ComputeCanonicalSubgraphAndList(int sortedIndex);

    template<typename T>
    std::vector<std::vector<T>> GetPowerSet(const std::vector<T>& elements); /// generate a power set of a given type

    bool IsEdgeInVertexSet(const UndirectedEdge& edge, const std::vector<int>& vertices); /// are the endpoints of edge in the given set of vertices

    void GenerateSubDiagrams(); /// main routine called by constructor

    void GenerateCanonicalSubDiagrams(); /// generate list of UNIQUE canonical subgraphs

    void GenerateEmbedListsForSubDiagrams(); /// generate VertexEmbedList objects for each subgraph

    GraphContainer CreateSubgraphContainerFromVertexMapAndEdges(const std::vector<UndirectedEdge>& edges, const std::vector<int>& map);

    int GetVertexSiteIndex(int vertexLabel) const; /// get lattice site index for a given vertex (through OriginalList)

    bool AreDisjoint(int index1, int index2); /// are two subgraphs, labeled by index1 and index2 disjoint

    void ComputeDisjointSets(); /// construct D_n's

    void AddToSortedSubdiagrams(const GraphContainer& g, int indexToUnsorted); /// add subgraph to the sorted set with index corresponding to the unsorted (VerticesMap)

    std::pair<int, int> IndexConversionSorted(int sortedIndex) const; /// sorted linear index 0,1,...,N_sub-1 mapped to (bondIndex,subIndex)

    std::vector<UndirectedEmbeddedEdge> ConvertUndirectedEdgesToUndirectedEmbeddedEdges(const std::vector<UndirectedEdge>& inputEdges); /// create embedded edges from edges

    int GetSortedIndexFromUnsortedIndex(int unsortedIndex) const;

public:
    SubDiagramGenerator(const GraphContainer& container, const VertexEmbedList& list, CubicLattice *lattice);

    void PrintSubDiagram(int index) const; /// print a given subdiagram

    void PrintSubDiagrams(); /// print subdiagrams (debugging purposes)

    static std::pair<std::vector<UndirectedEdge>, std::vector<int>> GetRelabeledEdgesAndVertexMap(const std::vector<UndirectedEdge>& edgeSet); /// relabel edges and provide map to original labels

    /**** public accessors ****/
    int GetNbrSubDiagrams() const { return this->NbrSubDiagrams; }

    int GetSizeDN(int n) const;

    int GetSortedLinearIndex(int nbrBonds, int graphIndex) const;

    int GetVertexMapIndexForSubDiagram(int sortedIndex) const; /// for a given sorted linear index, get the corresponding index for unsorted (VerticesMap and EmbedLists)

    GraphContainer GetSubDiagram(int sortedIndex) const;

    GraphContainer GetSubDiagram(int nbrBonds, int graphIndex) const;

    GraphContainer GetCanonicalSubDiagramContainer(int nbrBonds, int graphIndex) const;

    GraphContainer GetCanonicalSubDiagramContainer(int sortedIndex) const;

    VertexEmbedList GetCanonicalSubDiagramEmbedList(int nbrBonds, int graphIndex) const;

    VertexEmbedList GetCanonicalSubDiagramEmbedList(int sortedIndex) const;

    std::vector<int> GetVertexMap(int sortedIndex) const;

    std::vector<int> GetVertexMap(int nbrBonds, int graphIndex) const;

    std::vector<UndirectedEmbeddedEdge> GetEmbeddedEdgeSet(int sortedIndex) const;

    std::vector<UndirectedEmbeddedEdge> GetEmbeddedEdgeSet(int nbrBonds, int graphIndex) const;

    VertexEmbedList GetEmbedListCanonicalRelabel(int sortedIndex);

    VertexEmbedList GetEmbedListCanonicalRelabel(int nbrBonds, int graphIndex);

    VertexEmbedList GetEmbedList(int sortedIndex);

    VertexEmbedList GetEmbedList(int nbrBonds, int graphIndex) const;

    int GetSizeSubDiagrams(int nbrBonds) const;

};

#endif // SUBDIAGRAMGENERATOR_H
