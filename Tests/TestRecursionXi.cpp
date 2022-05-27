#include "GraphGeneratorNauty.h"
#include "CanonicalGraphManager.h"
#include "CubicLattice.h"
#include "VertexEmbedList.h"
#include "XiRecursion.h"
#include "XiRecursionRooted.h"

int main(int argc, char *argv[])
{
    int n = 4;
    int m = SETWORDSNEEDED(n);

    CanonicalGraphManager MyManager(n);
    std::cout << "AFTER_GRAPH_MANAGER!\n";

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainer(n, m, g, 2); /// container from densenauty
    TestContainer.SetRootedVertex(0, 0);
    TestContainer.SetRootedVertex(1, 1);
    std::cout << TestContainer;

    /// check that this is canonical!
    auto tempConContainer = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, std::vector<int>{0, 1});
    if (tempConContainer!=TestContainer)
        throw std::logic_error("ERROR: two colored canonical graphs not matching up!\n");

    CubicLattice MyLattice(100);

    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index2 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index3 = MyLattice.GetNearestNeighbor(index1, 2);
    auto index4 = MyLattice.GetNearestNeighbor(index3, 1);
    VertexEmbedList EmbedList(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor);
    EmbedList.AddFixedVertexEmbed(0, 1, index1);
    EmbedList.AddFixedVertexEmbed(1, 2, index2);
    EmbedList.AddVertexEmbed(3, index3);
    EmbedList.AddVertexEmbed(4, index4);

    std::cout << "PRINTING_PLAQ_XI\n";
    XiRecursionRooted<ZClusterPureGaugeArbEmbedding> RecursionObject(&MyManager, TestContainer, EmbedList, &MyLattice);
    std::cout << "PRINTING_PLAQ_FULL_XI_GINAC\n";
    std::cout << RecursionObject.GetFullXiGiNaC() << "\n";
    std::cout << "PRINTING_PLAQ_EXPANDED_XI_GINAC\n";
    std::cout << RecursionObject.GetExpandedXiGiNaC() << "\n";

    g6String = "A_"; /// graph
    tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer OneLinkContainer(2, m, g, 2); /// container from densenauty
    TestContainer.SetRootedVertex(0, 0);
    TestContainer.SetRootedVertex(1, 1);
    std::cout << OneLinkContainer;

    /// check that this is canonical!
    auto onelinkTempContainer = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, std::vector<int>{0, 1});
    if (onelinkTempContainer!=OneLinkContainer)
        throw std::logic_error("ERROR: two colored canonical graphs not matching up!\n");

    VertexEmbedList OneLinkEmbedList(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor);
    OneLinkEmbedList.AddFixedVertexEmbed(0, 1, index1);
    OneLinkEmbedList.AddFixedVertexEmbed(1, 2, index2);

    std::cout << "PRINTING_ONELINK_XI\n";
    XiRecursionRooted<ZClusterPureGaugeArbEmbedding> OneLinkRecursionObject(&MyManager, OneLinkContainer, OneLinkEmbedList, &MyLattice);
    std::cout << "PRINTING_ONELINK_FULL_XI_GINAC\n";
    std::cout << OneLinkRecursionObject.GetFullXiGiNaC() << "\n";
    std::cout << "PRINTING_ONELINK_EXPANDED_XI_GINAC\n";
    std::cout << OneLinkRecursionObject.GetExpandedXiGiNaC() << "\n";

    return 0;
}
