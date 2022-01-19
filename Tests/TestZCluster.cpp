#include <iostream>

#include "ZClusterPureGauge.h"
#include "ZClusterPureGaugeArbEmbedding.h"

int main()
{
    /******* test graph weights ********/
    //// "square": Cr (g6)
    int n = 4;
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// "square" graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    CubicLattice MyLattice(100);
    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex

    GraphContainer TestContainer(n, m, g); /// container from densenauty
    std::cout << TestContainer;

    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index3 = MyLattice.GetNearestNeighbor(index1, 2);
    auto index2 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index4 = MyLattice.GetNearestNeighbor(index2, 2);
    VertexEmbedList EmbedList(MaxInteractionLength::NearestNeighbor);
    EmbedList.AddVertexEmbed(1, index1);
    EmbedList.AddVertexEmbed(2, index2);
    EmbedList.AddVertexEmbed(3, index3);
    EmbedList.AddVertexEmbed(4, index4);

    ZClusterPureGauge MyCluster(&TestContainer, &EmbedList, &MyLattice);
    std::cout << "****SQUARE****\n";
    MyCluster.PrintZ();

    DYNFREE(g, g_sz); /// free graph

    /// "square with handle": DhS (g6)
    n = 5;
    m = SETWORDSNEEDED(n);

    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    g6String = "DhS"; /// "square with handle" graph
    tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainerHandle(n, m, g); /// container from densenauty
    std::cout << TestContainerHandle;

    index1 = MyLattice.GetSiteIndex(indicesCenter); /// embed with "long" handle (straight two-link coupling available!)
    index2 = MyLattice.GetNearestNeighbor(index1, 1);
    index3 = MyLattice.GetNearestNeighbor(index2, 1);
    index4 = MyLattice.GetNearestNeighbor(index3, 5);
    auto index5 = MyLattice.GetNearestNeighbor(index2, 5);
    VertexEmbedList EmbedListLongHanlde(MaxInteractionLength::NearestNeighbor);
    EmbedListLongHanlde.AddVertexEmbed(1, index1);
    EmbedListLongHanlde.AddVertexEmbed(2, index2);
    EmbedListLongHanlde.AddVertexEmbed(3, index3);
    EmbedListLongHanlde.AddVertexEmbed(4, index4);
    EmbedListLongHanlde.AddVertexEmbed(5, index5);

    ZClusterPureGauge MyClusterLongHandle(&TestContainerHandle, &EmbedListLongHanlde, &MyLattice);
    std::cout << "****LONG_HANDLE****\n";
    MyClusterLongHandle.PrintZ();

    index1 = MyLattice.GetSiteIndex(indicesCenter); /// embed with "bent" handle
    index2 = MyLattice.GetNearestNeighbor(index1, 3);
    index3 = MyLattice.GetNearestNeighbor(index2, 2);
    index4 = MyLattice.GetNearestNeighbor(index3, 1);
    index5 = MyLattice.GetNearestNeighbor(index2, 1);
    VertexEmbedList EmbedListBentHanlde(MaxInteractionLength::NearestNeighbor);
    EmbedListBentHanlde.AddVertexEmbed(1, index1);
    EmbedListBentHanlde.AddVertexEmbed(2, index2);
    EmbedListBentHanlde.AddVertexEmbed(3, index3);
    EmbedListBentHanlde.AddVertexEmbed(4, index4);
    EmbedListBentHanlde.AddVertexEmbed(5, index5);
    EmbedListBentHanlde.PrintList();

    ZClusterPureGauge MyClusterBentHandle(&TestContainerHandle, &EmbedListBentHanlde, &MyLattice);
    std::cout << "****BENT_HANDLE****\n";
    MyClusterBentHandle.PrintZ();

    index1 = MyLattice.GetSiteIndex(indicesCenter); /// embed with handle not being a nearest neighbor
    index2 = MyLattice.GetNextNearestNeighbor(index1, 1);
    index3 = MyLattice.GetNearestNeighbor(index2, 1);
    index4 = MyLattice.GetNearestNeighbor(index3, 2);
    index5 = MyLattice.GetNearestNeighbor(index2, 2);
    VertexEmbedList EmbedListNNNHandle(MaxInteractionLength::NextNearestNeighbor);
    EmbedListNNNHandle.AddVertexEmbed(1, index1);
    EmbedListNNNHandle.AddVertexEmbed(2, index2);
    EmbedListNNNHandle.AddVertexEmbed(3, index3);
    EmbedListNNNHandle.AddVertexEmbed(4, index4);
    EmbedListNNNHandle.AddVertexEmbed(5, index5);
    EmbedListNNNHandle.PrintList();

    ZClusterPureGaugeArbEmbedding MyClusterNNNHandle(&TestContainerHandle, &EmbedListNNNHandle, &MyLattice);
    std::cout << "****NNN_HANDLE****\n";
    MyClusterNNNHandle.PrintZ();

    DYNFREE(g, g_sz); /// free graph

    return 0;
}
