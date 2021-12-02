#include "VertexEmbedList.h"

/// constructor for correlator embedding (rooted graph)
/// takes in what the max embedding length will be and the correlator distance
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxEmbeddingLength, MaxInteractionLength correlatorDistance) :
    NbrChoicesForFirstBond(1),
    MaxLength(maxEmbeddingLength),
    TwoPointFunction(true),
    FixedVertices(2, VertexEmbed{-1,-1}),
    CorrelatorDistance(correlatorDistance)
{
#ifdef DEBUG
    std::cout << "Created a rooted graph for a correlator of size " << correlatorDistance << "!\n";
#endif
    switch (maxEmbeddingLength)
    {
    case MaxInteractionLength::NearestNeighbor:
        this->BondCounts.resize(1);
        break;
    case MaxInteractionLength::NextNearestNeighbor:
        this->BondCounts.resize(2);
        break;
    case MaxInteractionLength::ThirdNearestNeighbor:
        this->BondCounts.resize(3);
        break;
    case MaxInteractionLength::FourthNearestNeighbor:
        this->BondCounts.resize(4);
        break;
    default:
        throw std::invalid_argument("Invalid interaction length given to constructor of VertexEmbedList!\n");
        break;
    }
}

/// constructor for unrooted graph
/// takes in what the max interaction length will be
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxEmbeddingLength) :
    MaxLength(maxEmbeddingLength),
    TwoPointFunction(false),
    FixedVertices(2, VertexEmbed{-1,-1})
{
#ifdef DEBUG
    std::cout << "Created an unrooted VertexEmbedList!\n";
#endif
    switch (maxEmbeddingLength)
    {
    case MaxInteractionLength::NearestNeighbor:
        this->BondCounts.resize(1);
        break;
    case MaxInteractionLength::NextNearestNeighbor:
        this->BondCounts.resize(2);
        break;
    case MaxInteractionLength::ThirdNearestNeighbor:
        this->BondCounts.resize(3);
        break;
    case MaxInteractionLength::FourthNearestNeighbor:
        this->BondCounts.resize(4);
        break;
    default:
        throw std::invalid_argument("Invalid interaction length given to constructor of VertexEmbedList!\n");
        break;
    }
}

/// add a vertex to the list
void VertexEmbedList::AddVertexEmbed(const VertexEmbed& v)
{
#ifdef DEBUG
    for (auto it=this->List.begin(); it!=this->List.end(); ++it)
        if ((v.Number==it->Number) !=  (v.Index==it->Index))
            std::cout << "WARNING: AddVertexEmbed attempting to add VertexEmbed object with previously occupied site or previously used vertex label!\n";
#endif
    this->List.insert(v);
}

/// add a vertex to the list
void VertexEmbedList::AddVertexEmbed(int number, int index)
{
    return this->AddVertexEmbed(VertexEmbed{number,index});
}

/// set both fixed vertices (add entries to member variable List as well!)
void VertexEmbedList::AddFixedVerticesEmbed(const std::vector<VertexEmbed> &embed)
{
    if (embed.size() != 2)
        throw std::invalid_argument("AddFixedVerticesEmbed requires a vector of size 2!\n");
    if (this->FixedVertices[0].Index!=-1 || this->FixedVertices[1].Index!=-1)
        std::cerr << "WARNING: AddFixedVerticesEmbed fixed vertices already set!\n";
    if (!this->TwoPointFunction)
        std::cerr << "WARNING: AddFixedVerticesEmbed called but TwoPointFunction flag set to false!\n";
    this->FixedVertices = embed;
    this->List.insert(embed[0]);
    this->List.insert(embed[1]);
}

/// set fixed vertex labeled by fixedNbr (0 or 1) (add to member variable List as well!)
void VertexEmbedList::AddFixedVertexEmbed(int fixedNbr, const VertexEmbed& embed)
{
    if (fixedNbr != 0 && fixedNbr != 1)
        throw std::invalid_argument("AddFixedVertexEmbed requires fixedNbr to be 0 or 1!\n");
    if (this->FixedVertices[fixedNbr].Index!=-1)
        std::cerr << "WARNING: AddFixedVertexEmbed fixed vertices already set!\n";
    if (!this->TwoPointFunction)
        std::cerr << "WARNING: AddFixedVertexEmbed called but TwoPointFunction flag set to false!\n";
    this->FixedVertices[fixedNbr] = embed;
    this->List.insert(embed);
}

/// set fixed vertex labeled by fixedNbr (overload)
void VertexEmbedList::AddFixedVertexEmbed(int fixedNbr, int number, int index)
{
    this->AddFixedVertexEmbed(fixedNbr, VertexEmbed{number, index});
}

/// update the bond count
/// when a vertex is added, we need to check how many times to call this
/// will call atleast once as embedding algorithm will choose a neighbor (NN (0) or NNN (1) etc) of already placed vertex but new vertex could also be adjacent to other already placed vertices
void VertexEmbedList::IncrementBondCount(int dIndex)
{
    if (dIndex < 0 || dIndex >= this->BondCounts.size())
        throw std::invalid_argument("IncrementBondCount requires 0 <= dIndex < BondCounts.size()!\n");
    this->BondCounts[dIndex]++;
}

int VertexEmbedList::GetBondCount(int dIndex) const
{
    if (dIndex < 0 || dIndex >= this->BondCounts.size())
        throw std::invalid_argument("GetBondCount requires 0 <= dIndex < BondCounts.size()!\n");
    return this->BondCounts[dIndex];
}

/// check if list has repeated vertices
bool VertexEmbedList::HasRepeatedVertices() const
{
    for (auto it=this->List.begin(); it!=std::prev(this->List.end(),2); ++it)
        if (std::find_if(std::next(it), this->List.end(), [it](const VertexEmbed& v) { return (it->Number == v.Number); } ) != this->List.end())
            return true;
    return false;
}

/// check if list has repeated sites
bool VertexEmbedList::HasRepeatedSites() const
{
    for (auto it=this->List.begin(); it!=std::prev(this->List.end(),2); ++it)
        if (std::find_if(std::next(it), this->List.end(), [it](const VertexEmbed& v) { return (it->Index == v.Index); } ) != this->List.end())
            return true;
    return false;
}

/// accessor for fixed vertices
VertexEmbed VertexEmbedList::GetFixedVertex(int index) const
{
    if (!this->IsTwoPointFunction())
        throw std::logic_error("GetFixedVertex called for unrooted graph!\n");
    if (index !=0 && index != 1)
        throw std::invalid_argument("GetFixedVertex requires index to be 0 or 1!\n");
    return this->FixedVertices[index];
}

/// gets the color of the vertex: 0 for the first rooted vertex, 1 for the second rooted vertex, and 2 for all others (returns 2 for unrooted graphs)
/// @arg number: number of the vertex
int VertexEmbedList::GetVertexColor(int number) const
{
#ifdef DEBUG
    auto it = this->begin();
    for (; it!=this->end(); ++it)
        if (number==it->Number)
            break;
    if (it==this->end())
        std::cout << "ERR: GetVertexColor given a vertex number which is not in the list!\n";
#endif
    if (!this->IsTwoPointFunction())
        return 2;
    if (number==this->FixedVertices[0])
        return 0;
    if (number==this->FixedVertices[1])
        return 1;
    return 2;
}

std::vector<VertexEmbed> VertexEmbedList::GetSortedList() const
{
    return std::vector<VertexEmbed>(this->List.begin(), this->List.end());
}

/// get the site index (lattice) for a given vertex using it's label
int VertexEmbedList::GetVertexSiteIndex(int vertexLabel) const
{
    auto it = std::find_if(this->List.begin(), this->List.end(), [vertexLabel](const VertexEmbed& v) { return (v.Number==vertexLabel); });
    if (it==this->List.end()) /// could not find it
    {
        std::cerr << "ERROR: VertexEmbedList::GetVertexSiteIndex could not find the given vertex label!\n";
        return -1;
    }
    return it->Index;
}

/// equality operator for VertexEmbedList
bool operator==(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{

    if (lhs.IsTwoPointFunction()!=rhs.IsTwoPointFunction()) /// check if they are of the same type (rooted or unrooted)
        return false;

    if (lhs.GetSize()!=rhs.GetSize()) /// check that they are the same size
        return false;

    if (lhs.GetNbrBondTypes()!=rhs.GetNbrBondTypes()) /// check that they have the same number of bond types
        return false;

    for (int i=0; i<lhs.GetNbrBondTypes(); ++i) /// check bond counts
        if (lhs.GetBondCount(i)!=rhs.GetBondCount(i))
            return false;

    if (lhs.IsTwoPointFunction()) /// compare fixed vertices if we have a two-point function (order matters!)
        if ((lhs.GetFixedVertex(0) != rhs.GetFixedVertex(0)) || (lhs.GetFixedVertex(1) != rhs.GetFixedVertex(1)))
            return false;

    for (auto it1 = lhs.begin(), it2 = rhs.begin(); it1!=lhs.end() && it2!=rhs.end(); ++it1, ++it2) /// compare all of the embedded vertices (ordered)
    {
        if (*it1 != *it2)
            return false;
    }

    return true;
}

/// inequality operator
bool operator!=(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{
    return !(lhs == rhs);
}

/// less than operator which will be needed for std::set<VertexEmbedList>
/// lists with the smaller size come first
/// comparison is then made on this->List (lexigraphical order for vectors)
/// no reference is made to rooted or unrooted
bool operator<(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{
    if (lhs.GetSize()!=rhs.GetSize())
        return (lhs.GetSize()<rhs.GetSize());
    return (lhs.List < rhs.List);
}

/// output VertexEmbedList
std::ostream& operator<<(std::ostream& os, const VertexEmbedList& list)
{
    os << "********VertexEmbedList********\n";
    for (int i=0; i<list.GetNbrBondTypes(); ++i)
        os << "BondDegree " << i << " with count " << list.GetBondCount(i) << "\n";
    for (auto it=list.begin(); it!=list.end(); ++it)
        os << "Vertex " << *it << "\n";
    if (list.IsTwoPointFunction())
    {
        for (auto it=list.FixedVertices.begin(); it!=list.FixedVertices.end(); ++it)
            os << "ROOTED: " << *it << "\n";
    }
    return os;
}
