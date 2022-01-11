#include <iostream>
#include <chrono>

#include "GraphEmbedder.h"
#include "CubicLatticeCanonicalizor.h"
#include "GraphGeneratorNauty.h"

int main(int argc, char *argv[])
{
    GraphEmbedder MyEmbedder(argc, argv);
    auto resultCanonicalGraphsAndCounts = MyEmbedder.GetCanonicalGraphsAndCounts(); /// test this routine! (plaquette with handle!)
    for (int i=0; i<std::get<1>(resultCanonicalGraphsAndCounts).size(); ++i)
    {
        std::cout << "embedding " << i << " with counts " << std::get<2>(resultCanonicalGraphsAndCounts)[i] << "\n";
        std::cout << std::get<1>(resultCanonicalGraphsAndCounts)[i] << "\n";
    }
    auto resultPair = MyEmbedder.ContainerAndSampleCubicEmbeddingFromG6();

    CubicLattice MyLattice(100);
    CubicLatticeCanonicalizor MyCubicLatticeCanonicalizor1(&resultPair.first, &MyLattice, resultPair.second);
    CubicLatticeCanonicalizor MyCubicLatticeCanonicalizor2(&resultPair.first, &MyLattice, resultPair.second);

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

    if (resultOld==resultNew)
        std::cout << "Success! Explicit rotations and permutations agree!\n";
    else
        std::cout << "Failure! Explicit rotations and permutations do NOT agree!\n";

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

    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
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

    if (result1!=result2)
        std::cout << "RESULT1 and RESULT2 do not equal each other! Success!\n";
    else
        std::cout << "RESULT1 and RESULT2 equal each other! Failure!\n";

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

    /// canonicalize a rooted graph and its embeddings

    std::vector<int> rootedVertices1{0,2};
    auto rooted1 = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, rootedVertices1);
    std::cout << "TEST_CANONICAL_COLORED1:\n";
    rooted1.PrintM(); /***** C0---C1---C2 *****/

    std::vector<int> rootedVertices2{0,1};
    auto rooted2 = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, rootedVertices2);
    std::cout << "TEST_CANONICAL_COLORED2:\n";
    rooted2.PrintM(); /***** C0---C2---C1 *****/

    std::vector<int> rootedVertices3{2,0};
    auto rooted3 = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, rootedVertices3);
    std::cout << "TEST_CANONICAL_COLORED3:\n";
    rooted3.PrintM(); /***** C1---C0---C2 *****/

    index2 = MyLattice.GetNearestNeighbor(index1, 1);
    index3 = MyLattice.GetNearestNeighbor(index2, 1);
    VertexEmbedList EmbedListColor1Embed1(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// first embedding of first rooted grpah ("straight" embedding)
    EmbedListColor1Embed1.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor1Embed1.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor1Embed1.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor1Embed1(&rooted1, &MyLattice, EmbedListColor1Embed1);

    auto resultColor1Embed1 = TestCanonicalizorColor1Embed1.GetCanonical();
    std::cout << "RESULT_COLOR1_EMBED1: " << resultColor1Embed1 << "\n";

    index3 = MyLattice.GetNearestNeighbor(index2, 2);
    VertexEmbedList EmbedListColor1Embed2(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// second embedding of first rooted graph ("bent" in +y embedding)
    EmbedListColor1Embed2.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor1Embed2.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor1Embed2.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor1Embed2(&rooted1, &MyLattice, EmbedListColor1Embed2);

    auto resultColor1Embed2 = TestCanonicalizorColor1Embed2.GetCanonical();
    std::cout << "RESULT_COLOR1_EMBED2: " << resultColor1Embed2 << "\n";

    index3 = MyLattice.GetNearestNeighbor(index2, 5);
    VertexEmbedList EmbedListColor1Embed3(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// third embedding of first rooted graph ("bent" in -y embedding)
    EmbedListColor1Embed3.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor1Embed3.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor1Embed3.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor1Embed3(&rooted1, &MyLattice, EmbedListColor1Embed3);

    auto resultColor1Embed3 = TestCanonicalizorColor1Embed3.GetCanonical();
    std::cout << "RESULT_COLOR1_EMBED3: " << resultColor1Embed3 << "\n";

    if (resultColor1Embed2==resultColor1Embed3)
        std::cout << "SUCCESS! These rooted graph embeddings give same canonical graph!\n";
    else
        std::cout << "FAILURE! These rooted graph embeddings SHOULD give same canonilcal graph!\n";

    index2 = MyLattice.GetNearestNeighbor(index1, 4);
    index3 = MyLattice.GetNearestNeighbor(index1, 1);
    VertexEmbedList EmbedListColor3Embed1(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// first embedding of third rooted grpah ("straight" embedding)
    EmbedListColor3Embed1.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor3Embed1.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor3Embed1.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor3Embed1(&rooted3, &MyLattice, EmbedListColor3Embed1);

    auto resultColor3Embed1 = TestCanonicalizorColor3Embed1.GetCanonical();
    std::cout << "RESULT_COLOR3_EMBED1: " << resultColor3Embed1 << "\n";

    index2 = MyLattice.GetNearestNeighbor(index1, 4);
    index3 = MyLattice.GetNearestNeighbor(index1, 2);
    VertexEmbedList EmbedListColor3Embed2(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// second embedding of second rooted grpah ("bent" +y embedding)
    EmbedListColor3Embed2.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor3Embed2.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor3Embed2.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor3Embed2(&rooted3, &MyLattice, EmbedListColor3Embed2);

    auto resultColor3Embed2 = TestCanonicalizorColor3Embed2.GetCanonical();
    std::cout << "RESULT_COLOR3_EMBED2: " << resultColor3Embed2 << "\n";

    if (resultColor1Embed1 != resultColor1Embed2 && resultColor1Embed1 != resultColor3Embed1 && resultColor1Embed1 != resultColor3Embed2 && resultColor3Embed1 != resultColor3Embed2 && resultColor1Embed2 != resultColor3Embed1 && resultColor1Embed2 != resultColor3Embed2)
        std::cout << "SUCCESS! Embeddings of different rooted graphs not equal to each other!\n";
    else
        std::cout << "FAILURE! Embeddings of different rooted graphs SHOULD not equal to each other!\n";

    index3 = MyLattice.GetNearestNeighbor(index1, 1);
    index2 = MyLattice.GetNearestNeighbor(index3, 1);
    VertexEmbedList EmbedListColor2Embed1(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// first embedding of second rooted grpah ("straight" embedding)
    EmbedListColor2Embed1.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor2Embed1.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor2Embed1.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor2Embed1(&rooted2, &MyLattice, EmbedListColor2Embed1);

    auto resultColor2Embed1 = TestCanonicalizorColor2Embed1.GetCanonical();
    std::cout << "RESULT_COLOR2_EMBED1: " << resultColor2Embed1 << "\n";

    index3 = MyLattice.GetNearestNeighbor(index1, 1);
    index2 = MyLattice.GetNearestNeighbor(index3, 2);
    VertexEmbedList EmbedListColor2Embed2(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// second embedding of second rooted grpah ("bent" +y embedding)
    EmbedListColor2Embed2.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor2Embed2.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor2Embed2.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor2Embed2(&rooted2, &MyLattice, EmbedListColor2Embed2);

    auto resultColor2Embed2 = TestCanonicalizorColor2Embed2.GetCanonical();
    std::cout << "RESULT_COLOR2_EMBED2: " << resultColor2Embed2 << "\n";

    if (resultColor2Embed1!=resultColor2Embed2)
        std::cout << "SUCCESS! Embeddings of different rooted graphs not equal to each other!\n";
    else
        std::cout << "FAILURE! Embeddings of different rooted graphs SHOULD not equal to each other!\n";

    index3 = MyLattice.GetNearestNeighbor(index1, 1);
    index2 = MyLattice.GetNearestNeighbor(index3, 5);
    VertexEmbedList EmbedListColor2Embed3(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// second embedding of second rooted grpah ("bent" -y embedding)
    EmbedListColor2Embed3.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor2Embed3.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor2Embed3.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor2Embed3(&rooted2, &MyLattice, EmbedListColor2Embed3);

    auto resultColor2Embed3 = TestCanonicalizorColor2Embed3.GetCanonical();
    std::cout << "RESULT_COLOR2_EMBED3: " << resultColor2Embed3 << "\n";

    index3 = MyLattice.GetNearestNeighbor(index1, 1);
    index2 = MyLattice.GetNearestNeighbor(index3, 3);
    VertexEmbedList EmbedListColor2Embed4(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor); /// second embedding of second rooted grpah ("bent" +z embedding)
    EmbedListColor2Embed4.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    EmbedListColor2Embed4.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    EmbedListColor2Embed4.AddVertexEmbed(3, index3);

    CubicLatticeCanonicalizor TestCanonicalizorColor2Embed4(&rooted2, &MyLattice, EmbedListColor2Embed4);

    auto resultColor2Embed4 = TestCanonicalizorColor2Embed4.GetCanonical();
    std::cout << "RESULT_COLOR2_EMBED4: " << resultColor2Embed4 << "\n";

    if (resultColor2Embed2==resultColor2Embed3 && resultColor2Embed2==resultColor2Embed4)
        std::cout << "Success! These rooted graph embeddings are equivalent!\n";
    else
        std::cout << "Success! These embeddings SHOULD be equivalent!\n";

}
