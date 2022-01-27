#ifndef XIRECURSION_H
#define XIRECURSION_H

#include "CanonicalGraphManager.h"
#include "SubDiagramGenerator.h"

/// auxiliary class for the finite-cluster method
/// uniquely identifies terms which stand for log Z(g), which is the partition function on a given cluster g
class XiExpansionUnrootedTerm
{
private:
    std::pair<int,int> CanonicalCoordinates; /// coordinates for canonical container (unrooted): (bonds, index)

    std::vector<UndirectedEmbeddedEdge> EmbeddedVertices; /// embedded vertices (lattice indices only)

    VertexEmbedList FullyCanonicalEmbedList; /// canonical embed list (wrt to nauty AND cubic symmetries)

    int Coefficient; /// coefficient

public:
    XiExpansionUnrootedTerm(const std::pair<int,int>& coordinates, const std::vector<UndirectedEmbeddedEdge>& embeddedVertices, const VertexEmbedList& fullyCanonicalEmbedList, int coefficient) :
        CanonicalCoordinates(coordinates),
        EmbeddedVertices(embeddedVertices),
        FullyCanonicalEmbedList(fullyCanonicalEmbedList),
        Coefficient(coefficient)
    {}

    /***** public accessors *****/

    int GetCoefficient() const { return this->Coefficient; }

    std::pair<int,int> GetCanonicalCoordinates() const { return this->CanonicalCoordinates; }

    VertexEmbedList GetFullyCanonicalEmbedList() const { return this->FullyCanonicalEmbedList; }

    std::vector<UndirectedEmbeddedEdge> GetEmbeddedVertices() const { return this->EmbeddedVertices; }

    // coefficient \to coefficient + coefficientToAdd
    void CombineCoefficients(int coefficientToAdd) { this->Coefficient+=coefficientToAdd; }

    // user-defined operators
    friend std::ostream& operator<< (std::ostream& stream, const XiExpansionUnrootedTerm& can);
    friend bool operator==(const XiExpansionUnrootedTerm& lhs, const XiExpansionUnrootedTerm& rhs);

};

/// output (debugging purposes)
inline std::ostream& operator<<(std::ostream& stream, const XiExpansionUnrootedTerm& xi)
{
    stream << "XI_TERM:\n";
    stream << "(" << xi.CanonicalCoordinates.first << "," << xi.CanonicalCoordinates.second << ")\n";
    stream << xi.EmbeddedVertices << "\n";
    stream << xi.FullyCanonicalEmbedList;
    stream << "COEFFICIENT: " << xi.Coefficient << "\n";
    return stream;
}

/// equality and inequality operator for the terms in the xi expansion
inline bool operator==(const XiExpansionUnrootedTerm& lhs, const XiExpansionUnrootedTerm& rhs)
{
    /// first check the canonical coordinates (easy check so that we do not need to always compare edges)
    if (!((lhs.CanonicalCoordinates.first == rhs.CanonicalCoordinates.first) && (lhs.CanonicalCoordinates.second == rhs.CanonicalCoordinates.second)))
        return false;
    /// if they are equal then we compare embedded edges
    return (lhs.EmbeddedVertices==rhs.EmbeddedVertices);
}

inline bool operator!=(const XiExpansionUnrootedTerm& lhs, const XiExpansionUnrootedTerm& rhs)
{
    return !(lhs==rhs);
}

/// class which performs the Xi recursion given a graph and it's embedding (unrooted)
class XiRecursionUnrooted
{
private:
    CanonicalGraphManager *GraphManager; /// pointer to look-up table object with canonical graphs (rooted and unrooted)

    GraphContainer XiContainer; /// container for cluster

    VertexEmbedList XiEmbedList; /// embed list for cluster

    CubicLattice *Lattice; /// pointer to lattice object

    std::vector<XiExpansionUnrootedTerm> XiTerms; /// list of terms in the expansion for xi

    void ComputeXiTerms(GraphContainer container, VertexEmbedList embedList, int sign); /// recursive method for computing xi

    void AddXiTerm(const XiExpansionUnrootedTerm& newTerm); /// add term to list

public:
    XiRecursionUnrooted(CanonicalGraphManager* manager, const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice *lattice);

    VertexEmbedList RestoreOriginalLabels(const VertexEmbedList& embedList); /// restore original vertex labels

    /**** public accessors ****/

    int GetNbrXiTerms() const { return this->XiTerms.size(); }

    XiExpansionUnrootedTerm GetXiTerm(int index) const;

};

#endif // XIRECURSION_H
