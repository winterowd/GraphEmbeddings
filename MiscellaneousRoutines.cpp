#include "MiscellaneousRoutines.h"

int OppositeDir(const int& dir)
{
    return (dir+2)%6+1;
}

/// save this in miscellaneous cpp file!
void PegsInHoles(const unsigned int n, const unsigned int m, bool identicalPegs, bool verbose)
{
    if (m>n)
        throw std::invalid_argument("PegsInHoles: M must be less than or equal to N!");

    std::vector<int> s(m, -1);
    std::vector<bool> occ(n, false);

    int count = 0;
    int sk;
    int k=0;
    while (1)
    {
        if (s[k] < 0 && k > 0 && identicalPegs)
            sk = s[k-1]+1;
        else
            sk = s[k]+1; /// next location of kth peg
        s[k] = sk; /// set kth peg

        if (sk>=n) /// check if we went over the number of holes
        {
            /// reset kth peg
            s[k] = -1; /// distinguishable pegs
            k--; /// go back to previous peg (k-1)
            if (k<0) /// we are at the end if the first peg has reached past the last hole
                break;
            if (verbose)
                std::cout << "Went over number of holes! Marking hole " << s[k] << " occupied by peg " << k << " as unoccupied!\n";
            occ[s[k]] = false; /// mark as unoccupied as we intend to move (k-1) peg to new location
        }
        else
        {
            if (!occ[sk]) /// check if occupied
            {
                if (verbose)
                    std::cout << "Setting peg " << k << " at " << sk << "!\n";
                occ[sk] = true; /// set occupied if it is
                k++; /// go to next peg
            }
            else
            {
                if (verbose)
                    std::cout << "Tried setting peg " << k << " at " << sk << " which was occupied! " << occ[sk] << "\n";
                continue; /// if it is occupied, move kth peg again
            }

            if (k>=m) /// have we successfully placed all the pegs?
            {
                count++;
                std::cout << "FOUND A CONFIG! " << k << " " << count << std::endl;
                for (int i=0; i<m; ++i)
                    std::cout << s[i] << " ";
                std::cout << "\n";
                for (int i=0; i<n; ++i)
                    std::cout << occ[i] << " ";
                std::cout << "\n";

                k--; /// go back (we are over the size of the array s)

                if (verbose)
                    std::cout << "Marking hole " << s[k] << " occupied by peg " << k << " as unoccupied!\n";

                occ[s[k]] = false; /// mark as unoccupied as we intend to move kth peg to new location
            }

        }

    }

}

bool DoesNotDoubleBack(const std::vector<int> &path)
{
    for (auto i=0; i<path.size()-1; ++i)
    {
        if (path[i]==OppositeDir(path[i+1]))
            return false;
    }
    return true;
}

void GenerateAllPermutationsWithRepeats(std::vector<std::vector<int>> &lists, const std::string& s, std::vector<int>& pos, int n, const int& size)
{
    if (n == size)
    {
        std::vector<int> temp;
        for (int i=0; i<n; ++i)
            temp.push_back(std::stoi(std::string(1, s[pos[i]])));
        if (DoesNotDoubleBack(temp))
            lists.push_back(temp);
        return;
    }
    for (int i=0; i<s.size(); ++i)
    {
        pos[n] = i;
        GenerateAllPermutationsWithRepeats(lists, s, pos, n+1, size);
    }
}

bool operator==(const Site& lhs, const Site& rhs)
{
    return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

bool operator!=(const Site& lhs, const Site& rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<< (std::ostream& stream, const Site& site)
{
    stream << site.x << " " << site.y << " " << site.z;
    return stream;
}

Site MoveDir(const Site& start, int dir, int L)
{
    auto result = Site(start);
    switch (dir) {
    case 1:
        result.x = (start.x+1)%L;
        break;
    case 2:
        result.y = (start.y+1)%L;
        break;
    case 3:
        result.z = (start.z+1)%L;
        break;
    case 4:
        result.x = (result.x+L-1)%L;
        break;
    case 5:
        result.y = (result.y+L-1)%L;
        break;
    case 6:
        result.z = (result.z+L-1)%L;
        break;
     default:
        throw std::invalid_argument("Invalid nearest neighbor index for CubicLattice!");
    }
    return result;
}

/// save this in miscellaneous cpp file!
/// is the bridge self-avoiding and does it end on the correct site
bool IsBridgeValid(const Site& start, const Site& end, const std::vector<int>& bridgeDirs, std::vector<Site>& Bridge, int L, bool verbose)
{
    Bridge.push_back(start);
    for (unsigned int i=0; i<bridgeDirs.size(); ++i)
    {
        auto newSite = MoveDir(Bridge[i], bridgeDirs[i], L);
        if (std::find(Bridge.begin(), Bridge.end(), newSite) != Bridge.end())
        {
            if (!(start == end && i==(bridgeDirs.size()-1))) /// can be a single loop!
            {
                if (verbose)
                {
                    std::cout << "******BRIDGE IS NOT SELF-AVOIDING!******" << std::endl;
                    Bridge.push_back(newSite);
                    for (auto j=0; j<Bridge.size(); ++j)
                        std::cout << Bridge[j] << std::endl;
                    std::cout << "****************************************" << std::endl;
                }
                return false;
            }
        }
        Bridge.push_back(newSite);
    }
    if (Bridge.back()!=end)
        return false;
    return true;
}

void TestLexicographicalOrderingVertexEmbedList()
{
    std::set<VertexEmbedList> setOfVertexEmbedLists;

    /***** create a "star" on the square lattice *****/
    SquareLattice lattice(100); /// lattice object
    VertexEmbedList firstVertexList(MaxInteractionLength::NearestNeighbor); /// first VertexEmbedList

    std::vector<unsigned int> initialSiteIndices(lattice.GetDim(), lattice.GetN()/2);
    auto initialIndex = lattice.GetSiteIndex(initialSiteIndices);
    firstVertexList.AddVertexEmbed(VertexEmbed{2,initialIndex});

    auto tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 1);
    firstVertexList.AddVertexEmbed(VertexEmbed{4,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 2);
    firstVertexList.AddVertexEmbed(VertexEmbed{1,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 3);
    firstVertexList.AddVertexEmbed(VertexEmbed{5,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 4);
    firstVertexList.AddVertexEmbed(VertexEmbed{3,tempIndex});

    setOfVertexEmbedLists.insert(firstVertexList);

    VertexEmbedList secondVertexList(MaxInteractionLength::NearestNeighbor); /// second VertexEmbedList

    initialIndex = lattice.GetSiteIndex(initialSiteIndices);
    secondVertexList.AddVertexEmbed(VertexEmbed{2,initialIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 1);
    secondVertexList.AddVertexEmbed(VertexEmbed{4,tempIndex-22}); /// smaller site index than in first list (lexigraphical order should sort this!)

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 2);
    secondVertexList.AddVertexEmbed(VertexEmbed{1,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 3);
    secondVertexList.AddVertexEmbed(VertexEmbed{5,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 4);
    secondVertexList.AddVertexEmbed(VertexEmbed{3,tempIndex});

    setOfVertexEmbedLists.insert(secondVertexList);

    VertexEmbedList thirdVertexList(MaxInteractionLength::NearestNeighbor); /// third VertexEmbedList (3 vertices)

    initialIndex = lattice.GetSiteIndex(initialSiteIndices);
    thirdVertexList.AddVertexEmbed(VertexEmbed{2,initialIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 1);
    thirdVertexList.AddVertexEmbed(VertexEmbed{4,tempIndex-22});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 2);
    thirdVertexList.AddVertexEmbed(VertexEmbed{1,tempIndex});

    setOfVertexEmbedLists.insert(thirdVertexList);

    VertexEmbedList fourthVertexList(MaxInteractionLength::NearestNeighbor); /// fourth VertexEmbedList

    initialIndex = lattice.GetSiteIndex(initialSiteIndices);
    fourthVertexList.AddVertexEmbed(VertexEmbed{2,initialIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 1);
    fourthVertexList.AddVertexEmbed(VertexEmbed{4,tempIndex-22}); /// smaller site index than in first list (lexigraphical order should sort this!)

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 4);
    fourthVertexList.AddVertexEmbed(VertexEmbed{3,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 3);
    fourthVertexList.AddVertexEmbed(VertexEmbed{5,tempIndex});

    tempIndex = lattice.GetNextNearestNeighbor(initialIndex, 2);
    fourthVertexList.AddVertexEmbed(VertexEmbed{1,tempIndex});

    setOfVertexEmbedLists.insert(fourthVertexList);

    std::cout << "VERTEX_EMBED_LIST has size: " << setOfVertexEmbedLists.size() << "\n";
    std::cout << "TestSetofVertexEmbedLists initial...\n";
    for (auto it=setOfVertexEmbedLists.begin(); it!=setOfVertexEmbedLists.end(); ++it)
        std::cout << " " << *it << "\n";

    int vertexCount = 5; /// get rid of VertexEmbedList objects that do not have 5 elements...
    auto it = setOfVertexEmbedLists.begin();
    while (it!=setOfVertexEmbedLists.end())
    {
        if (it->GetSize()!=vertexCount)
            it = setOfVertexEmbedLists.erase(it);
        else
            ++it;
    }

    std::cout << "TestSetofVertexEmbedLists after erasing...\n";
    for (auto it=setOfVertexEmbedLists.begin(); it!=setOfVertexEmbedLists.end(); ++it)
        std::cout << " " << *it << "\n";

}
