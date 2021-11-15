#include <iostream>
#include <chrono>

#include "GraphEmbedder.h"
#include "CubicLatticeCanonicalizor.h"

int main(int argc, char *argv[])
{
    GraphEmbedder MyEmbedder(argc, argv);
    auto resultPair = MyEmbedder.ContainerAndSampleCubicEmbeddingFromG6();
    CubicLattice MyLattice(100);
    CubicLatticeCanonicalizor MyCubicLatticeCanonicalizor1(&resultPair.first, &MyLattice, resultPair.second);
    CubicLatticeCanonicalizor MyCubicLatticeCanonicalizor2(&resultPair.first, &MyLattice, resultPair.second);
    std::vector<double> coords{1.1,2.2,3.3};
    //MyCubicLatticeCanonicalizor.TestPermutationsAsOrthogonalGroup(coords);
    auto start = std::chrono::high_resolution_clock::now();
    auto resultOld = MyCubicLatticeCanonicalizor1.GetCanonicalOld();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "DURATION_OLD: " << duration.count() << "\n";
    std::cout << "DEBUG_RESULTOLD: " << resultOld;
    start = std::chrono::high_resolution_clock::now();
    auto resultNew = MyCubicLatticeCanonicalizor2.GetCanonical();
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "DURATION_NEW: " << duration.count() << "\n";
    std::cout << "DEBUG_RESULTNEW: " << resultNew;

    /********** test **********/

    int n = 3;
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "BW"; /// graph: 1---3---2
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer TestContainer(n, m, g); /// container from densenauty

    std::vector<unsigned int> indices1(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indices1); /// index for vertex 1
    auto index3 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index2 = MyLattice.GetNearestNeighbor(index3, 1);

    /***** test embeddings ******/

    VertexEmbedList EmbedList1(MaxInteractionLength::NearestNeighbor); /// first embedding ("straight" embedding)
    EmbedList1.AddVertexEmbed(1, index1);
    EmbedList1.AddVertexEmbed(3, index3);
    EmbedList1.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor1(&TestContainer, &MyLattice, EmbedList1);

    auto result1 = TestCanonicalizor1.GetCanonical();
    std::cout << "RESULT1: " << result1 << "\n";

    index2 = MyLattice.GetNearestNeighbor(index3, 2); // change vertex 2 only

    VertexEmbedList EmbedList2(MaxInteractionLength::NearestNeighbor); ///  ("bent upwards in y" embedding)
    EmbedList2.AddVertexEmbed(1, index1);
    EmbedList2.AddVertexEmbed(3, index3);
    EmbedList2.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor2(&TestContainer, &MyLattice, EmbedList2);

    auto result2 = TestCanonicalizor2.GetCanonical();
    std::cout << "RESULT2: " << result2 << "\n";

    index2 = MyLattice.GetNearestNeighbor(index3, 5); /// change vertex 2 only

    VertexEmbedList EmbedList3(MaxInteractionLength::NearestNeighbor); /// ("bent downwards in y" embedding)
    EmbedList3.AddVertexEmbed(1, index1);
    EmbedList3.AddVertexEmbed(3, index3);
    EmbedList3.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor3(&TestContainer, &MyLattice, EmbedList3);

    auto result3 = TestCanonicalizor3.GetCanonical();
    std::cout << "RESULT3: " << result3 << "\n";

    index2 = MyLattice.GetNearestNeighbor(index3, 3); /// change vertex 2 only

    VertexEmbedList EmbedList4(MaxInteractionLength::NearestNeighbor); ///  ("bent upwards in z" embedding)
    EmbedList4.AddVertexEmbed(1, index1);
    EmbedList4.AddVertexEmbed(3, index3);
    EmbedList4.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor4(&TestContainer, &MyLattice, EmbedList4);

    auto result4 = TestCanonicalizor4.GetCanonical();
    std::cout << "RESULT4: " << result4 << "\n";

    index2 = MyLattice.GetNearestNeighbor(index3, 6); /// change vertex 2 only

    VertexEmbedList EmbedList5(MaxInteractionLength::NearestNeighbor); ///  ("bent upwards in z" embedding)
    EmbedList5.AddVertexEmbed(1, index1);
    EmbedList5.AddVertexEmbed(3, index3);
    EmbedList5.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor5(&TestContainer, &MyLattice, EmbedList5);

    auto result5 = TestCanonicalizor5.GetCanonical();
    std::cout << "RESULT5: " << result5 << "\n";

    index3 = MyLattice.GetNearestNeighbor(index1, 4);
    index2 = MyLattice.GetNearestNeighbor(index3, 6);

    VertexEmbedList EmbedList6(MaxInteractionLength::NearestNeighbor); /// some other
    EmbedList6.AddVertexEmbed(1, index1);
    EmbedList6.AddVertexEmbed(3, index3);
    EmbedList6.AddVertexEmbed(2, index2);

    CubicLatticeCanonicalizor TestCanonicalizor6(&TestContainer, &MyLattice, EmbedList6);

    auto result6 = TestCanonicalizor6.GetCanonical();
    std::cout << "RESULT6: " << result6 << "\n";
    result6 = TestCanonicalizor6.GetCanonical();
    std::cout << "RESULT6: " << result6 << "\n";

    if (result2==result3 && result2==result4 && result2==result5 && result2==result6)
        std::cout << "RESULT2, RESULT3, RESULT4, RESULT5, RESULT6 are identical! Success!\n";
    else
        std::cout << "RESULT2, RESULT3, RESULT4, RESULT5, RESULT6 are NOT identical! Failure!\n";


    DYNFREE(g, g_sz); /// free graph

}
