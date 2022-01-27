#include "GraphGeneratorNauty.h"
#include "CanonicalGraphManager.h"
#include "CubicLattice.h"
#include "VertexEmbedList.h"
#include "XiRecursion.h"

int main(int argc, char *argv[])
{
    /*int n = 3;
    int m = SETWORDSNEEDED(n);

    CanonicalGraphManager MyManager(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "BW"; /// graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainer(n, m, g); /// container from densenauty
    std::cout << "DEBUG_THREE_LINK:\n";
    std::cout << TestContainer;

    CubicLattice MyLattice(100);

    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index2 = MyLattice.GetThirdNearestNeighbor(index1, 1);
    auto index3 = MyLattice.GetNextNearestNeighbor(index1, 1);
    VertexEmbedList EmbedList(MaxInteractionLength::NearestNeighbor);
    EmbedList.AddVertexEmbed(1, index1);
    EmbedList.AddVertexEmbed(2, index2);
    EmbedList.AddVertexEmbed(3, index3);

    XiRecursion RecursionObject(&MyManager, TestContainer, EmbedList, &MyLattice);*/

    int n = 4;
    int m = SETWORDSNEEDED(n);

    CanonicalGraphManager MyManager(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainer(n, m, g); /// container from densenauty
    //std::cout << "DEBUG_THREE_LINK:\n";
    std::cout << TestContainer;

    CubicLattice MyLattice(100);

    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index2 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index3 = MyLattice.GetNearestNeighbor(index1, 2);
    auto index4 = MyLattice.GetNearestNeighbor(index3, 1);
    VertexEmbedList EmbedList(MaxInteractionLength::NearestNeighbor);
    EmbedList.AddVertexEmbed(1, index1);
    EmbedList.AddVertexEmbed(2, index2);
    EmbedList.AddVertexEmbed(3, index3);
    EmbedList.AddVertexEmbed(4, index4);

    XiRecursionUnrooted RecursionObject(&MyManager, TestContainer, EmbedList, &MyLattice);

    return 0;
}
