#ifndef SUBDIAGRAMGENERATOR_H
#define SUBDIAGRAMGENERATOR_H

#include <algorithm>

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "CubicLattice.h"
#include "CubicLatticeCanonicalizor.h"

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

/// helper class in order to distinguish canonical vertex lists with different number of bonds
/// TODO: move this to CubicLatticeCanonicalizor and change return types etc so that this is returned
class CanonicalSubDiagram
{
private:
    int N; /// number of vertices

    int L; /// number of bonds

    VertexEmbedList CanonicalList; /// canonical vertex embed list (output from CubicLatticeCanonicalizor)
public:
    CanonicalSubDiagram(int l, const VertexEmbedList& list) : N(list.GetSize()), L(l), CanonicalList(list) {}

    int GetN() const { return this->N; }

    int GetL() const { return this->L; }

    friend bool operator==(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs);
    friend std::ostream& operator<< (std::ostream& stream, const CanonicalSubDiagram& can);
};

/// equality operator (need to first compare number of bonds and vertices!)
inline bool operator==(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs)
{
    if (lhs.GetN()!=rhs.GetN())
        return false;
    if (lhs.GetL()!=rhs.GetL())
        return false;
    return (lhs.CanonicalList==rhs.CanonicalList);
}

inline bool operator!=(const CanonicalSubDiagram& lhs, const CanonicalSubDiagram& rhs)
{
   return !(lhs==rhs);
}

inline std::ostream& operator<< (std::ostream& stream, const CanonicalSubDiagram& can)
{
    std::cout << "CANONICAL_GRAPH: L " << can.L << " N " << can.N << "\n";
    std::cout << can.CanonicalList;
    return stream;
}

/// generate all subgraphs of a given graph (embedded on a lattice i.e. cubic)
class SubDiagramGenerator
{
private:
    /***** private variables *****/
    GraphContainer *OriginalContainer; /// pointer to original graph container

    VertexEmbedList *OriginalList; /// pointer to original vertex embed list

    CubicLattice *MyCubicLattice; /// pointer to cubic lattice object

    std::vector<std::vector<std::pair<int, GraphContainer>>> SortedSubDiagramsWithMap; /// subdiagrams sorted by number of bonds together with (UNSORTED) index matching the vertex map

    std::vector<std::vector<int>> VerticesMap; /// map of vertices of subgraphs (order differs with SortedSubDiagrams)

    std::vector<VertexEmbedList> EmbedLists; ///  embed list for the subgraphs (same order as VerticesMap)

    std::vector<int> SubgraphToCanonicalMap;

    /// TODO: add EmbedLists that are canonical with respect to NAUTY and the octohedral group! (cubic symmetry)
    std::vector<CanonicalSubDiagram> CanonicalEmbedLists;

    std::vector<SetOfSets<int>> DisjointSets; /// D_n: set of subsets represented by the index corresponding to ordering of subgraphs in SortedSubDiagramsWithMap
    int NbrSubDiagrams; /// total number of connected subdiagrams

    /**** private methods ****/

    //// canonicalize with respect to NAUTY and cubic symmetries
    CanonicalSubDiagram ComputeCanonicalSubgraph(int sortedIndex);

    template<typename T>
    std::vector<std::vector<T>> GetPowerSet(const std::vector<T>& elements); /// generate a power set of a given type

    bool IsEdgeInVertexSet(const UndirectedEdge& edge, const std::vector<int>& vertices); /// are the endpoints of edge in the given set of vertices

    void GenerateSubDiagrams(); /// main routine called by constructor

    void GenerateCanonicalSubDiagrams(); /// generate list of UNIQUE canonical subgraphs

    void GenerateEmbedListsForSubDiagrams(); /// generate VertexEmbedList objects for each subgraph

    int GetVertexSiteIndex(int vertexLabel) const; /// get lattice site index for a given vertex (through OriginalList)

    bool AreDisjoint(int index1, int index2); /// are two subgraphs, labeled by index1 and index2 disjoint

    void ComputeDisjointSets(); /// construct D_n's

    void AddToSortedSubdiagrams(const GraphContainer& g, int indexToUnsorted); /// add subgraph to the sorted set with index corresponding to the unsorted (VerticesMap)

    std::pair<int, int> IndexConversionSorted(int sortedIndex) const; /// sorted linear index 0,1,...,N_sub-1 mapped to (bondIndex,subIndex)

public:
    SubDiagramGenerator(GraphContainer *container, VertexEmbedList *list, CubicLattice *lattice);

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

    std::vector<int> GetVertexMap(int sortedIndex) const;

    std::vector<int> GetVertexMap(int nbrBonds, int graphIndex) const;

    VertexEmbedList GetEmbedList(int sortedIndex) const;

    VertexEmbedList GetEmbedList(int nbrBonds, int graphIndex) const;

    int GetSizeSubDiagrams(int nbrBonds) const;

};

#endif // SUBDIAGRAMGENERATOR_H
