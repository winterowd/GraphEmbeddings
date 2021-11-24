#ifndef VERTEXEMBEDLIST_H
#define VERTEXEMBEDLIST_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "AbstractLattice.h"

/// data structure for an embedded vertex
struct VertexEmbed {
    int Number; /// vertex number
    int Index; /// lattice index
};

inline std::ostream& operator<<(std::ostream& os, const VertexEmbed& v)
{
    os << v.Number << " at site " << v.Index;
    return os;
}

/// equality and inequality based on Index and Number
inline bool operator==(const VertexEmbed& lhs, const VertexEmbed& rhs)
{
    return (lhs.Index == rhs.Index && lhs.Number == rhs.Number);
}

inline bool operator!=(const VertexEmbed& lhs, const VertexEmbed& rhs)
{
    return !(lhs==rhs);
}

/// equality and inequality based on only Index
inline bool operator==(const VertexEmbed& lhs, const int& rhs)
{
    return (lhs.Index == rhs);
}

inline bool operator!=(const VertexEmbed& lhs, const int& rhs)
{
    return !(lhs == rhs);
}

/// and repeat for the Index appearing on the left hand side...
inline bool operator==(const int& lhs,  const VertexEmbed& rhs)
{
    return (lhs == rhs.Index);
}

inline bool operator!=(const int& lhs, const VertexEmbed& rhs)
{
    return !(lhs == rhs);
}

/// used for sorting VertexEmbed objects!
/// compare vertex numbers, if equal then compare indices
inline bool operator<(const VertexEmbed& lhs,  const VertexEmbed& rhs)
{
    if (lhs.Number != rhs.Number)
        return (lhs.Number < rhs.Number);
    else
        return (lhs.Index < rhs.Index);
}

/// class to hold lists for embedding
/// when including next nearest neighbors and so on need to keep track of number of number of types of links used
/// getters and setters mostly. What other types of routines? any other data fields?
class VertexEmbedList
{
private:
    std::set<VertexEmbed> List; /// data for embedded graphs!

    std::vector<int> BondCounts; /// only need to keep track of total counts and not which edges are connected by NN or NNN

    int NbrChoicesForFirstBond; /// remember number of choices for first bond for unrooted graphs or 1 for rooted graphs

    MaxInteractionLength MaxLength; /// save this

    bool TwoPointFunction; /// flag for two-point function (rooted graph)

     /// NOTE: when making canonical form, first vertex corresponds to color 1 and second vertex to color 2 (equality operator will distinguish them)
    std::vector<VertexEmbed> FixedVertices; /// for two-point functions (rooted graph)

    MaxInteractionLength CorrelatorDistance; /// keep track of which type of correlator

public:

    VertexEmbedList(MaxInteractionLength maxLength); /// constructor for unrooted graph

    VertexEmbedList(MaxInteractionLength maxLength, MaxInteractionLength correlatorDistance); /// constructor for correlator (rooted graph)

    VertexEmbedList(const VertexEmbedList& list) = default;

    MaxInteractionLength GetMaxLength() const { return this->MaxLength; }

    void AddVertexEmbed(const VertexEmbed& v);

    void AddVertexEmbed(int number, int index);

    void AddFixedVerticesEmbed(const std::vector<VertexEmbed>& embed);

    void AddFixedVertexEmbed(int fixedNbr, const VertexEmbed& embed);

    void AddFixedVertexEmbed(int fixedNbr, int number, int index);

    void IncrementBondCount(int dIndex);

    int GetBondCount(int dIndex) const;

    VertexEmbed GetFixedVertex(int index) const;

    int GetVertexColor(int number) const;

    /**** accessors ****/
    int GetSize() const { return this->List.size(); }

    int GetNbrBondTypes() const { return this->BondCounts.size(); }

    bool IsTwoPointFunction() const { return this->TwoPointFunction; }

    MaxInteractionLength GetCorrelatorDistance() const { return this->CorrelatorDistance; }

    int GetCorrelatorDistanceAsIndex() const { return static_cast<int>(this->CorrelatorDistance); }

    bool HasRepeatedVertices() const; /// debugging routine

    bool HasRepeatedSites() const; /// debugging routine

    void SetNbrChoicesForFirstBond(int nbr) { this->NbrChoicesForFirstBond = nbr; }

    int GetNbrChoicesForFirstBond() const { return this->NbrChoicesForFirstBond; }

    std::vector<VertexEmbed> GetSortedList() const; /// return sorted list for comparisons

    void PrintList() const { std::cout << "VertexEmbedList::PrintList():\n"; for (auto it=this->List.begin(); it!=this->List.end(); ++it) std::cout << *it << "\n"; }

    /// define iterator types
    using iterator = std::set<VertexEmbed>::iterator;
    using const_iterator = std::set<VertexEmbed>::const_iterator;

    /// begin and end functions
    iterator begin() { return this->List.begin(); }
    const_iterator begin() const { return this->List.begin(); }
    iterator end() { return this->List.end(); }
    const_iterator end() const { return this->List.end(); }

    friend std::ostream& operator<<(std::ostream& os, const VertexEmbedList& list); /// output for debugging reasons
    friend bool operator<(const VertexEmbedList& lhs, const VertexEmbedList& rhs);

};

bool operator==(const VertexEmbedList& lhs, const VertexEmbedList& rhs); /// comparison operator
bool operator!=(const VertexEmbedList& lhs, const VertexEmbedList& rhs); /// comparison operator

#endif // VERTEXEMBEDLIST_H
