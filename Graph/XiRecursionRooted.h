#ifndef XIRECURSIONROOTED_H
#define XIRECURSIONROOTED_H

#include "CanonicalGraphManager.h"
#include "SubDiagramGenerator.h"
#include "TwoPointCorrelator.h"

/// auxiliary class for the finite-cluster method
/// uniquely identifies terms which stand for log Z(g), which is the partition function on a given cluster g
class XiExpansionRootedTerm
{
private:
    std::pair<int,int> CanonicalCoordinates; /// coordinates for canonical container (rooted): (bonds, index)

    std::vector<UndirectedEmbeddedEdge> EmbeddedVertices; /// embedded vertices (lattice indices only)

    VertexEmbedList FullyCanonicalEmbedList; /// canonical embed list (wrt to nauty AND cubic symmetries)

    VertexEmbedList EmbedListOriginal; /// embed list where rooted vertices reside at ORIGINAL location

    int Coefficient; /// coefficient

public:
    XiExpansionRootedTerm(const std::pair<int,int>& coordinates, const std::vector<UndirectedEmbeddedEdge>& embeddedVertices, const VertexEmbedList& fullyCanonicalEmbedList, const VertexEmbedList& originalEmbedList, int coefficient) :
        CanonicalCoordinates(coordinates),
        EmbeddedVertices(embeddedVertices),
        FullyCanonicalEmbedList(fullyCanonicalEmbedList),
        EmbedListOriginal(originalEmbedList),
        Coefficient(coefficient)
    {
        if (FullyCanonicalEmbedList.GetNbrSetRootedVertices()!=2)
            throw std::invalid_argument("Error: XiExpansionRootedTerm requires fullyCanonicalEmbedList to have two rooted vertices!\n");

        if (EmbedListOriginal.GetNbrSetRootedVertices()!=2)
            throw std::invalid_argument("Error: XiExpansionRootedTerm requires originalEmbedList to have two rooted vertices!\n");
    }

    /***** public accessors *****/

    int GetCoefficient() const { return this->Coefficient; }

    std::pair<int,int> GetCanonicalCoordinates() const { return this->CanonicalCoordinates; }

    VertexEmbedList GetFullyCanonicalEmbedList() const { return this->FullyCanonicalEmbedList; }

    std::vector<UndirectedEmbeddedEdge> GetEmbeddedVertices() const { return this->EmbeddedVertices; }

    // coefficient \to coefficient + coefficientToAdd
    void CombineCoefficients(int coefficientToAdd) { this->Coefficient+=coefficientToAdd; }

    // user-defined operators
    friend std::ostream& operator<< (std::ostream& stream, const XiExpansionRootedTerm& can);
    friend bool operator==(const XiExpansionRootedTerm& lhs, const XiExpansionRootedTerm& rhs);

};

/// output (debugging purposes)
inline std::ostream& operator<<(std::ostream& stream, const XiExpansionRootedTerm& xi)
{
    stream << "XI_TERM:\n";
    stream << "(" << xi.CanonicalCoordinates.first << "," << xi.CanonicalCoordinates.second << ")\n";
    stream << xi.EmbeddedVertices << "\n";
    stream << xi.FullyCanonicalEmbedList;
    stream << xi.EmbedListOriginal;
    stream << "COEFFICIENT: " << xi.Coefficient << "\n";
    return stream;
}

/// equality and inequality operator for the terms in the xi expansion
inline bool operator==(const XiExpansionRootedTerm& lhs, const XiExpansionRootedTerm& rhs)
{
    /// first check the canonical coordinates (easy check so that we do not need to always compare edges)
    if (!((lhs.CanonicalCoordinates.first == rhs.CanonicalCoordinates.first) && (lhs.CanonicalCoordinates.second == rhs.CanonicalCoordinates.second)))
        return false;
    /// TODO: do we also compare FullyCanonicalEmbedList?
    /// if they are equal then we compare embedded edges
    return (lhs.EmbeddedVertices==rhs.EmbeddedVertices);
}

inline bool operator!=(const XiExpansionRootedTerm& lhs, const XiExpansionRootedTerm& rhs)
{
    return !(lhs==rhs);
}

class XiRecursionRooted
{
private:
    CanonicalGraphManager *GraphManager; /// pointer to look-up table object with canonical graphs (rooted and unrooted)

    GraphContainer XiContainer; /// container for cluster

    VertexEmbedList XiEmbedList; /// embed list for cluster

    CubicLattice *Lattice; /// pointer to lattice object

    int EmbeddingNumber; /// embedding number

    std::vector<XiExpansionRootedTerm> XiTerms; /// list of terms in the expansion for xi

    std::vector<TwoPointCorrelator> CorrelatorTerms;

    void ComputeXiTerms(GraphContainer container, VertexEmbedList embedList, int sign); /// recursive method for computing xi

    void AddXiTerm(const XiExpansionRootedTerm& newTerm); /// add term to list

public:
    XiRecursionRooted(CanonicalGraphManager* manager, const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice *lattice, int embeddingNumber=1);

    /**** public accessors ****/

    int GetNbrXiTerms() const { return this->XiTerms.size(); }

    XiExpansionRootedTerm GetXiTerm(int index) const;
};

#endif // XIRECURSIONROOTED_H
