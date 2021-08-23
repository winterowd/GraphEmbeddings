#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <unordered_set>
#include <numeric>

extern "C" {
#include "nausparse.h"
}

#include "GenericConnectedGraph.h"
#include "GraphGenerator.h"
#include "GraphGeneratorNauty.h"
#include "GraphContainer.h"
#include "GraphEmbedder.h"
#include "SquareLattice.h"

/// DEBUG: global variables
int L = 10; /// lattice size
std::vector<std::vector<int>> BridgeDirections;

struct Site {
    int x, y, z;
};

bool operator==(const Site& lhs, const Site& rhs)
{
    return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

bool operator!=(const Site& lhs, const Site& rhs)
{
    return !(lhs == rhs);
}

/// save this in miscellaneous cpp file!
void PegsInHoles(const unsigned int n, const unsigned int m, bool identicalPegs, bool verbose=false)
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

/// DEBUG: global variables
Site MyStart{4,4,4}; /// starting site
Site MyEnd{4,4,4}; /// ending site

std::ostream& operator<< (std::ostream& stream, const Site& site)
{
    stream << site.x << " " << site.y << " " << site.z;
    return stream;
}

Site MoveDir(const Site& start, int dir)
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

int OppositeDir(const int& dir) { return (dir+2)%6+1; }

bool DoesNotDoubleBack(const std::vector<int> &path)
{
    for (auto i=0; i<path.size()-1; ++i)
    {
        if (path[i]==OppositeDir(path[i+1]))
            return false;
    }
    return true;
}

/// save this in miscellaneous cpp file!
/// is the bridge self-avoiding and does it end on the correct site
bool IsBridgeValid(const Site& start, const Site& end, const std::vector<int>& bridgeDirs, std::vector<Site>& Bridge, bool verbose=false)
{
    Bridge.push_back(start);
    for (unsigned int i=0; i<bridgeDirs.size(); ++i)
    {
        auto newSite = MoveDir(Bridge[i], bridgeDirs[i]);
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

/// save this in miscellaneous cpp file!
void GenerateAllPermutationsWithRepeats(const std::string& s, std::vector<int>& pos, int n, const int& size)
{
    if (n == size)
    {
        std::vector<int> temp;
        for (int i=0; i<n; ++i)
            temp.push_back(std::stoi(std::string(1, s[pos[i]])));
        if (DoesNotDoubleBack(temp))
            BridgeDirections.push_back(temp);
        return;
    }
    for (int i=0; i<s.size(); ++i)
    {
        pos[n] = i;
        GenerateAllPermutationsWithRepeats(s, pos, n+1, size);
    }
}

int CodeFromBool(const std::vector<bool>& vec)
{
    int result = 0;
    for (int i=0; i<vec.size(); ++i)
        result += vec[i] << i;
    return result;
}

void TestBreakFromFor(int& row, int& col)
{
    for (int i=0; i<5; ++i)
    {
        for (int j=i+1; j<5; ++j)
        {
            if ((i+j)%5==0)
            {
                row = i;
                col = j;
                break;
            }
        }
    }
    std::cout << "exiting at " << row << " " << col << "\n";
}

void TestVectorOfPairs(const std::vector<std::vector<bool>>& m, unsigned int n)
{
    std::vector<std::pair<int, int>> vertices;
    for (int i=0; i<n-1; i++) /// assume n <= m.size()
    {
        for (int j=i+1; j<n; j++) /// assume square
        {
            if (m[i][j])
                vertices.push_back(std::make_pair(i+1,j+1));
        }
    }
    for (int i=0; i<vertices.size(); ++i)
        std::cout << "vertex " << i+1 << ": " << vertices[i].first << " " << vertices[i].second << "\n";
}

void TestSparseNauty(const std::vector<std::vector<bool>>& adj)
{
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    SG_DECL(sg);
    SG_DECL(csg);
    statsblk stats;

    static DEFAULTOPTIONS_SPARSEGRAPH(options);

    options.getcanon = TRUE;

    int n = adj.size();
    int m = SETWORDSNEEDED(n);

    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

    auto nde = 0;
    for (int i=0; i<adj.size(); ++i)
        for (int j=i+1; i<adj.size(); ++i)
            if (adj[i][j])
                nde++;

    SG_ALLOC(sg, n, nde, "malloc");
    SG_ALLOC(csg, n, nde, "malloc");

    sg.nv = n;
    sg.nde = nde;

    sg.v[0] = 0;
    for (int i=0; i<adj.size(); ++i) /// set vertices and their degrees
    {
        //sg.v[i] = i;
        auto count = 0;
        for (int j=0; j<adj.size(); ++j)
        {
            if (adj[i][j])
                count++;
        }
        sg.d[i] = count;
        if (i>0 && count>0)
            sg.v[i] = sg.v[i-1]+sg.d[i];
    }

    for (int i=0; i<adj.size(); ++i)
    {
        if (sg.d[i]>0)
        {
            auto count = 0;
            for (int j=0; j<adj.size(); ++j)
            {
                if (adj[i][j])
                {
                    sg.e[sg.v[i]+count] = j;
                    count++;
                }

            }
            if (count != sg.d[i])
            {
                std::cerr << "ERROR! count does not equal vertex " << i << " order! " << count << " " << sg.d[i] << std::endl;
                return;
            }
        }
    }

    sparsenauty(&sg, lab, ptn, orbits, &options, &stats, &csg);

    if (aresame_sg(&csg,&sg))
        std::cout << "Graph is canonical!" << std::endl;
    else
        std::cout << "Graph is not canonical!" << std::endl;

    for (int i=0; i<adj.size(); ++i)
    {
        std::cout << "Vertex " << i << " connected to: ";
        for (int j=0; j<csg.d[i]; ++j)
        {
            std::cout << csg.e[csg.v[i]+j] << " ";
        }
        std::cout << "\n";
    }

}

void TestDenseNauty()
{
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    DYNALLSTAT(graph, g, g_sz); /// graph
    DYNALLSTAT(graph, cg, cg_sz); /// canonical graph
    statsblk stats;

    static DEFAULTOPTIONS_GRAPH(options); /// options

    options.getcanon = TRUE; /// set optoins

    int n = 4;
    int m = SETWORDSNEEDED(n); /// array of m setwords sufficients to hold n bits

    std::cout << "n: " << n << " m: " << m << std::endl;

    nauty_check(WORDSIZE, n, m, NAUTYVERSIONID); /// check if everything is ok

    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");

    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// graph
    DYNALLOC2(graph, cg, cg_sz, n, m, "malloc"); /// canonical graph

    EMPTYGRAPH(g, n, m); /// clear graph

    ADDONEEDGE(g, 0, 1, m);
    ADDONEEDGE(g, 0, 2, m);
    ADDONEEDGE(g, 0, 3, m);
    ADDONEEDGE(g, 1, 2, m);

    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);

    set *gj, *cgj;

    for (int j = 1; j < n; ++j)
    {
        gj = GRAPHROW(g,j,m);
        cgj = GRAPHROW(cg,j,m);
        for (int i = 0; i < j; ++i)
        {
            if (ISELEMENT(gj,i))
                std::cout << "G: Vertex " << j << " adjacent to vertex " << i << "\n";
            if (ISELEMENT(cgj,i))
                std::cout << "CG: Vertex " << j << " adjacent to vertex " << i << "\n";
        }
    }

}

void separate(std::string text, std::string s)
{
    std::cout << "BEFORE: " << text << " " << text.size() << std::endl;
    for (auto x : s)
        std::replace(text.begin(), text.end(), x, ' ');
    std::cout << "AFTER: " << text << " " << text.size() << std::endl;

    std::vector<std::string> tokens;
    std::size_t start = 0, end = 0;
    while ((end = text.find(' ', start)) != std::string::npos)
    {
        if (end != start)
        {
            tokens.push_back(text.substr(start, end-start));
            std::cout << "FOUND TOKEN! " << text.substr(start, end-start) << " " << end-start << "\n";
        }
        start = end + 1;
    }
    if (end != start)
    {
        tokens.push_back(text.substr(start));
        std::cout << "FOUND TOKEN! " << text.substr(start) << "\n";
    }
}

void TestContainer()
{
    /***** TEST FOR NEW CONTAINER ********/

    DYNALLSTAT(graph, g, g_sz); /// graph

    int n = 4;
    int m = SETWORDSNEEDED(n); /// array of m setwords sufficients to hold n bits

    std::cout << "n: " << n << " m: " << m << std::endl;

    nauty_check(WORDSIZE, n, m, NAUTYVERSIONID); /// check if everything is ok

    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// graph

    ADDONEEDGE(g, 0, 1, m);
    ADDONEEDGE(g, 0, 2, m);
    ADDONEEDGE(g, 0, 3, m);
    ADDONEEDGE(g, 1, 2, m);

    GraphContainer MyContainer(n,m,g);
}

void CheckErase()
{
    int size = 3;
    std::vector<std::vector<int>> vec;
    vec.push_back(std::vector<int>{-1,2});
    vec.push_back(std::vector<int>{3,4,10});
    vec.push_back(std::vector<int>{11,25});
    vec.erase(std::remove_if( vec.begin(), vec.end(), [&size](const std::vector<int> &p) { return p.size()<size; }), vec.end());
    for (int i=0; i<vec.size(); ++i)
        for (int j=0; j<vec[i].size(); ++j)
            std::cout << "AFTER_REMOVE: list " << i << " " << " element " << j << " " << vec[i][j] << "\n";
}

void CheckEqualityVertexEmbed()
{
    std::vector<VertexEmbed> vertextList;
    vertextList.push_back(VertexEmbed{1,25}); /// vertex 1 at lattice site index 25
    vertextList.push_back(VertexEmbed{2,2}); /// vertex 2 at lattice site index 2

    VertexEmbed tempVertex{1,25};

    if (std::find(vertextList.begin(), vertextList.end(), tempVertex) != vertextList.end())
        std::cout << "SUCCESS! FOUND ENTIRE VERTEX IN VERTEX LIST!\n";

    if (tempVertex != vertextList[0])
        std::cout << "SOMETHING WRONG WITH FULL VERTEX EQUALITY!\n";

    if (std::find(vertextList.begin(), vertextList.end(), 25) != vertextList.end())
        std::cout << "SUCCESS! FOUND LATTICE INDEX IN VERTEX LIST!\n";

}

void TestBreakFor()
{
    int i;
    for (i=0; i<5; ++i)
    {
        if (i==3)
            break;
    }
    std::cout << "TestBreakFor: " << i << "\n";
}

void TestAddingToVector()
{
    std::vector<int> vec{1,3,5,7};

    auto currentSize = vec.size();
    for (int i=0; i<currentSize; ++i) /// vector increases in size: what happens?
        vec.push_back(vec[i]*2);

    /// remove add elements
    vec.erase(std::remove_if(vec.begin(), vec.end(), [](const int& x) { return (x%2==1); }), vec.end());

    std::cout << "RESULT: ";
    for (int i=0; i<vec.size(); ++i)
            std::cout << vec[i] << " ";
    std::cout << "\n";
}

void TestGetRemainingVertices()
{
    std::unordered_set<int> mySetVertices{4,3,2,1}; /// precomputed in constructor of GraphEmbedder

    std::vector<VertexEmbed> vertextEmbedList; /// precomputed and passed in by embedding function
    vertextEmbedList.push_back(VertexEmbed{1,25}); /// vertex 1 at lattice site index 25
    vertextEmbedList.push_back(VertexEmbed{2,2}); /// vertex 2 at lattice site index 2
    vertextEmbedList.push_back(VertexEmbed{3,69}); /// vertex 3 at lattice site index 69

    /**** HERE STARTS THE FUNCTION THAT WOULD GO INTO GraphEmbedder ****/
    std::unordered_set<int> usedVertices;
    std::for_each(vertextEmbedList.begin(), vertextEmbedList.end(), [&usedVertices](const VertexEmbed& v) { usedVertices.insert(v.Number); });

    std::unordered_set<int> remainingVertices;

    for (std::unordered_set<int>::iterator it=mySetVertices.begin(); it!=mySetVertices.end(); ++it)
    {
        auto search = usedVertices.find(*it);
        if (search == usedVertices.end())
            remainingVertices.insert(*it);
    }

    for (auto x: remainingVertices)
        std::cout << "RemainingVertices " << x << "\n";

}

void TestFindIdenticalGraph()
{
    std::vector<std::vector<VertexEmbed>> lists;

    std::vector<VertexEmbed> first{VertexEmbed{1,20}, VertexEmbed{2,24}};
    lists.push_back(first);

    std::vector<VertexEmbed> second{VertexEmbed{1,20}, VertexEmbed{2,24}, VertexEmbed{3,40}, VertexEmbed{4,41}};
    lists.push_back(second);

    std::vector<VertexEmbed> third{VertexEmbed{1,20}, VertexEmbed{2,24}, VertexEmbed{3,40}, VertexEmbed{4,66}};
    lists.push_back(third);

    std::vector<VertexEmbed> toAdd{VertexEmbed{1,20}, VertexEmbed{2,24}, VertexEmbed{3,40}, VertexEmbed{4,41}};

    bool toBedAdded = true;

    for (int i=0; i<lists.size(); ++i)
    {
        if (lists[i].size() == toAdd.size())
        {
            auto ve = std::next(toAdd.begin(),2);
            for ( ; ve!=toAdd.end(); ++ve)
            {
                std::cout << "Checking element in toAdd against list " << i << "\n";
                if (std::find(std::next(lists[i].begin(), 2), lists[i].end(), *ve) == lists[i].end())
                {
                    std::cout << "toAdd contains VertexEmbed object with vertex " << ve->Number << " at site " << ve->Index << " which is not in list " << i << "\n";
                    break;
                }
                else
                {
                    std::cout << "toAdd contains VertexEmbed object with vertex " << ve->Number << " at site " << ve->Index << " which is also in list " << i << "\n";
                }
            }
            if (ve==toAdd.end())
            {
                toBedAdded = false;
                std::cout << "toAdd identical to list " << i << "\n";
                break;
            }
        }
    }
    if (toBedAdded)
        std::cout << "List can be added! No duplicates!\n";
}

/// pass by address and is modified
void TestReference2(int *x)
{
    *x = -99;
}

/// pass by reference where the reference is passed by address
void TestReference1(int &x)
{
    std::cout << "BEFORE: " << x << "\n";
    TestReference2(&x);
    std::cout << "AFTER: " << x << "\n";
}

class AbstractBase
{
private:
    int val;

protected:
    int GetVal() const { return val; }

public:
    AbstractBase(int x) :val(x) {}

    virtual ~AbstractBase() {};
};

class Descendant : public AbstractBase
{
public:
    Descendant(int x) : AbstractBase(x) {}

    ~Descendant() {}

    int GetValFromBase() const { return this->GetVal(); }

};

enum SomeEnum
{
    First,
    Second,
    Third
};

/// for pegs in holes
unsigned int N = 4;
unsigned int M = 2;

/// TODO: two executables!
/// 1) generates
/// 2) one that reads in (need to make parameters for that! BOOST!)

void TestRepeats()
{
    std::vector<int> tempList{1,23,23,23,35,35};
    std::vector<int>::iterator it;
    for (int i=0; i<(tempList.size()-1); ++i)
    {
        it = std::find(std::next(tempList.begin(), i+1), tempList.end(), tempList[i]);
        if (it != tempList.end())
            std::cout << "FOUND_REPEAT! element " << i << " = " << tempList[i] << " same as element " << *it << "\n";
    }
    //std::cout << "NO_REPEATS!\n";
    //return false;
}

void TestVertexEmbedListFunctions()
{
    /// do some tests here!
    VertexEmbedList tempList(MaxInteractionLength::NearestNeighbor);

    tempList.AddVertexEmbed(VertexEmbed{1,24});
    tempList.AddVertexEmbed(VertexEmbed{2,45});
    tempList.AddVertexEmbed(VertexEmbed{1,55});
    tempList.AddVertexEmbed(VertexEmbed{3,45});

    if (tempList.HasRepeatedVertices())
        std::cout << "ERROR: Repeated vertices found!\n";

    if (tempList.HasRepeatedSites())
        std::cout << "ERROR: Repeated sites found!\n";

    auto temp = std::find(tempList.begin(), tempList.end(), 55);
    if (temp!=tempList.end())
        std::cout << "ERROR: Lattice site not free! " << *temp << "\n";

    std::vector<VertexEmbedList> lists;
    lists.push_back(tempList);

    if (GraphEmbedder::IsDuplicate(lists, tempList))
        std::cout << "ERROR: GraphEmbedder found a duplicate!\n";

    VertexEmbedList tempList2(MaxInteractionLength::NearestNeighbor);
    tempList2.AddVertexEmbed(VertexEmbed{1,24});
    tempList2.AddVertexEmbed(VertexEmbed{2,45});
    tempList2.AddVertexEmbed(VertexEmbed{4,55});
    tempList2.AddVertexEmbed(VertexEmbed{3,45});

    if (!GraphEmbedder::IsDuplicate(lists, tempList2))
        std::cout << "ERROR: GraphEmbedder did not find a duplicate!\n";

    lists.push_back(tempList2);

    VertexEmbedList tempList3(MaxInteractionLength::NearestNeighbor);
    tempList3.AddVertexEmbed(VertexEmbed{1,24});
    tempList3.AddVertexEmbed(VertexEmbed{2,45});
    tempList3.AddVertexEmbed(VertexEmbed{4,55});
    tempList3.IncrementBondCount(0);
    tempList3.IncrementBondCount(0);
    tempList3.IncrementBondCount(0);

    lists.push_back(tempList3);

    std::cout << "BEFORE_ERASING: " << lists.size() << "\n";
    for (auto v: lists)
        std::cout << v << "\n";
    GraphEmbedder::TestEraseWrongSizes(lists, 4);
    std::cout << "AFTER_ERASING: " << lists.size() << "\n";
    for (auto v: lists)
        std::cout << v << "\n";

    /// TODO: test copy constructor of VertexEmbedList!
    std::cout << "TEST_COPY_CTOR_1: " << tempList3 << "\n";
    VertexEmbedList tempList4(tempList3);
    tempList3.IncrementBondCount(0);
    std::cout << "TEST_COPY_CTOR_2: " << tempList4 << "\n";
    std::cout << "TEST_COPY_CTOR_3: " << tempList3 << "\n";

}

std::vector<std::vector<int>> combos;

void Combinations(const std::vector<int>& arr /* input array */, std::vector<int>& data /* stores current combo */, int index)
{
    if (index == data.size()) /// created a combo?
    {
        if (std::accumulate(data.begin(), data.end(), 0)==arr.back()) /// combo satisfies \sum_i data_i = arr.back()?
        {
            combos.push_back(data);
        }
        return;
    }

    for (int i=0; i<arr.size(); ++i)
    {
        data[index] = arr[i];
        Combinations(arr, data, index+1);
    }
}

int main(int argc, char *argv[])
{

    //TestRepeats();
    TestVertexEmbedListFunctions();
    /*std::vector<int> arr{0,1,2,3,4,5};
    std::vector<int> data(3,-1);
    Combinations(arr, data, 0);

    for (int i=0; i<combos.size(); ++i)
    {
        std::cout << "COMBO " << i << "\n";
        for (auto x: combos[i])
            std::cout << x << " ";
        std::cout << "\n";
    }*/

    return 0;

    //std::string sep = " ,/";
    //std::string tokens = "first_token second_token,/ third_token/fourth_token";
    //separate(tokens, sep);

    //Descendant myDescendant(3);
    //std::cout << "GetValFromBase(): " << myDescendant.GetValFromBase() << "\n";
    //SomeEnum x;
    //std::cout << "AFTER DECLARING ENUM: " << x << "\n";
    //x = SomeEnum::Second;
    //std::cout << "AFTER ASSIGNING IT A VALUE: " << x << "\n";
    //return 0;
    //GraphGeneratorNauty MyGenerator(argc, argv);
    //MyGenerator.Generate();

    //CheckErase();
    //CheckEqualityVertexEmbed();
    //TestGetRemainingVertices();
    //TestBreakFor();
    //TestAddingToVector();
    //TestFindIdenticalGraph();
    //int x = 22;
    //TestReference1(x);
    //std::cout << "IN_MAIN: " << x << "\n";
    //return 0;

    //SquareLattice MyLattice(100);

    GraphEmbedder MyEmbedder(argc, argv);
    //return 0;
    MyEmbedder.Embed();

    return 0;
}
