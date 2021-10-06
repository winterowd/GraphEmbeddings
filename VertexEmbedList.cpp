#include "VertexEmbedList.h"

/// constructor for correlator embedding (rooted graph)
/// takes in what the max embedding length will be and the correlator distance
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxEmbeddingLength, MaxInteractionLength correlatorDistance) :
    NbrChoicesForFirstBond(1),
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
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxLength) :
    TwoPointFunction(false),
    FixedVertices(2, VertexEmbed{-1,-1})
{
#ifdef DEBUG
    std::cout << "Created an unrooted graph!\n";
#endif
    switch (maxLength)
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
    this->List.push_back(v);
}

/// add a vertex to the list
void VertexEmbedList::AddVertexEmbed(int number, int index)
{
    this->List.push_back(VertexEmbed{number,index});
}

void VertexEmbedList::AddFixedVerticesEmbed(const std::vector<VertexEmbed> &embed)
{
    if (embed.size() != 2)
        throw std::invalid_argument("AddFixedVerticesEmbed requires a vector of size 2!\n");
    if (this->FixedVertices[0].Index!=-1 || this->FixedVertices[1].Index!=-1)
        std::cerr << "WARNING: AddFixedVerticesEmbed fixed vertices already set!\n";
    if (!this->TwoPointFunction)
        std::cerr << "WARNING: AddFixedVerticesEmbed called but TwoPointFunction flag set to false!\n";
    this->FixedVertices = embed;
    this->List.push_back(embed[0]);
    this->List.push_back(embed[1]);
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

/// return a vertex in list by reference
VertexEmbed VertexEmbedList::GetVertexEmbed(int index) const
{
    if (index < 0 || index >= this->List.size())
        throw std::invalid_argument("GetVertexEmbed requires 0 <= index < List.size()!\n");
    return this->List[index];
}

/// check if list has repeated vertices
bool VertexEmbedList::HasRepeatedVertices()
{
    for (int i=0; i<(this->List.size()-1); ++i)
        if (std::find_if(std::next(this->List.begin(), i+1), this->List.end(), [this, i](const VertexEmbed& v) { return (this->List[i].Number == v.Number); } ) != this->List.end())
            return true;
    return false;
}

/// check if list has repeated sites
bool VertexEmbedList::HasRepeatedSites()
{
    for (int i=0; i<(this->List.size()-1); ++i)
        if (std::find(std::next(this->List.begin(), i+1), this->List.end(), this->List[i].Index) != this->List.end())
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

    for (auto it = rhs.begin(); it!=rhs.end(); ++it) /// compare all of the vertices
        if (std::find(lhs.begin(), lhs.end(), *it) == lhs.end())
            return false;

    return true;
}

/// inequality operator
bool operator!=(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{
    return !(lhs == rhs);
}

/// output VertexEmbedList
std::ostream& operator<<(std::ostream& os, const VertexEmbedList& list)
{
    os << "********VertexEmbedList********\n";
    for (int i=0; i<list.GetNbrBondTypes(); ++i)
        os << "BondDegree " << i << " with count " << list.GetBondCount(i) << "\n";
    for (int i=0; i<list.GetSize(); ++i)
        os << "Vertex " << list.GetVertexEmbed(i) << "\n";
    return os;
}
