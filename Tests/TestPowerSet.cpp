#include <set>
#include <vector>
#include <iostream>

#include "SubDiagramGenerator.h"
#include "CubicLattice.h"

template<typename T>
std::vector<std::vector<T>> GetPowerSet(const std::vector<T>& elements)
{
    if (elements.empty())
    {
        //std::cout << "RETURNING THE EMPTY SET!\n";
        return std::vector<std::vector<T>>(1, std::vector<T>());
    }
    auto allOfThem = GetPowerSet(std::vector<T>(elements.begin()+1,elements.end()));
    T elem = elements[0];
    //std::cout << "APPENDING " << elem << " to all sets in current power set!\n";
    const int n = allOfThem.size();
    //std::cout << "n: " << n << "\n";
    for (int i=0; i<n; ++i)
    {
        const std::vector<T>& s = allOfThem[i];
        allOfThem.push_back(s);
        allOfThem.back().push_back(elem);
    }
    //std::cout << "RETURNING NON-EMPTY-SET!\n";
    return allOfThem;
}

template<typename T>
std::set<T> DummyFunc(const std::set<T>& in)
{
    return in;
}

int main()
{
    std::vector<int> dummySet{1,2,3};

    auto result = GetPowerSet(dummySet);
    for (auto setIt=result.begin(); setIt!=result.end(); ++setIt)
    {
        std::cout << "[ ";
        for (auto elemIt=setIt->begin(); elemIt!=setIt->end(); ++elemIt)
            std::cout << *elemIt << ", ";
        std::cout << "]\n";
    }

    /******* test graph weights ********/
    int n = 4;
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// "square" graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainer(n, m, g); /// container from densenauty
    std::cout << TestContainer;

    CubicLattice MyLattice(100);

    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index3 = MyLattice.GetNearestNeighbor(index1, 2);
    auto index2 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index4 = MyLattice.GetNearestNeighbor(index2, 2);
    VertexEmbedList EmbedList(MaxInteractionLength::NearestNeighbor);
    EmbedList.AddVertexEmbed(1, index1);
    EmbedList.AddVertexEmbed(2, index2);
    EmbedList.AddVertexEmbed(3, index3);
    EmbedList.AddVertexEmbed(4, index4);

    SubDiagramGenerator MySubDiagrams(&TestContainer, &EmbedList, &MyLattice);
    std::cout << "PRINTING_SUBDIAGRAMS!\n";
    MySubDiagrams.PrintSubDiagrams();

    DYNFREE(g, g_sz); /// free graph

    return 0;
}
