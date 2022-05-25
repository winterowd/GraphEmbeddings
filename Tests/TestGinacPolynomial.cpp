#include "AuxiliaryRoutinesForGINAC.h"
#include "MyLambdaPolynomial.h"
#include "PureGaugeweight.h"
#include "PureGaugeWeightOld.h"
#include "ZClusterPureGaugeArbEmbedding.h"
#include "StaticQuarkWeight.h"
#include "TwoPointCorrelator.h"
#include "XiRecursionRooted.h"
#include "StaticQuarkWeight.h"

int main()
{
    std::vector<std::pair<double, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients1;
    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients2;
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder1{1,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder2{3,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder3{1,1,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder4{2,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder5{0,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder6{0,1,0,0};

    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1., tempOrder1});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{3., tempOrder2});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1./6., tempOrder3});

    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder5});
    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder4});
    //coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1,6), tempOrder3});

    auto result1 = AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients<double>(coefficients1, 6);
    auto result2 = AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients<GiNaC::numeric>(coefficients2, 6);

    //std::cout << "RESULT1: " << result1 << "\n";
    //std::cout << "RESULT2: " << result2 << "\n";

    // add two polynomials and see if result is correct
    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients3;
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder1});
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(3), tempOrder2});
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1,6), tempOrder3});
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(10), tempOrder4});
    coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder6}); /// counted twice! should simplify to 20 \lambda^2_1
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder5});
    //coefficients3.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder6});

    auto expandedTest = AuxiliaryRoutinesForGinac::GetLambdaExpandedPureGaugeRationalFunction(coefficients3, coefficients2, 8);
    std::cout << "EXPANDED_TEST: " << expandedTest << " " << GiNaC::is_a<GiNaC::pseries>(expandedTest) << " " << GiNaC::is_a<GiNaC::ex>(expandedTest) << "\n";

    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder7{0,0,0,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder8{0,0,0,0,1,0};
    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder9{0,0,0,0,0,1};

    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>> coefficients4;
    coefficients4.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(1), tempOrder7});
    coefficients4.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(2), tempOrder8});
    coefficients4.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(3), tempOrder9});
    auto someTempSingleFlavor = GiNaC::numeric(1)/AuxiliaryRoutinesForGinac::SingleFlavorPolynomialFromCoefficients(coefficients4, 8);

    std::cout << "Test_someTempSingleFlavor: " << someTempSingleFlavor << "\n";

    auto someResult = AuxiliaryRoutinesForGinac::GetSingleFlavorH1AndHBar1ExpandedFunction(someTempSingleFlavor, 2, 2);

    std::cout << "Test_someResult: " << someResult << "\n";

    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder10{0,0,0,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder11{1,0,0,0,1,0};
    std::array<int, MaxInteractionLength::NbrInteractions+2> tempOrder12{1,1,0,0,0,1};

    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>> coefficients5;
    coefficients5.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(1), tempOrder10});
    coefficients5.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(2), tempOrder11});
    coefficients5.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(3), tempOrder12});
    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>> coefficients6;
    coefficients6.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions+2>>{GiNaC::numeric(1), tempOrder10});
    auto resultFullExpansion = AuxiliaryRoutinesForGinac::GetFullyExpandedSingleFlavorRationalFunction(coefficients6, coefficients5, 4, 2, 2);

    std::cout << "Test_resultFullExpansion: " << resultFullExpansion << "\n";

    //return 0;

    /*** SOME_TESTS_FOR_GINAC::NUMERIC ***/
    int n = 4;
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// "square" graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer WeightContainer(n, m, g, 2); /// container from densenauty
    /// vertices 1 and 2 are rooted
    WeightContainer.SetRootedVertex(0,0);
    WeightContainer.SetRootedVertex(1,1);
    std::cout << WeightContainer;
    std::vector<ExternalPolyakovLoop> rootedVertices{ExternalPolyakovLoop{1,true}, ExternalPolyakovLoop{2,false}};
    std::vector<ExternalPolyakovLoopOld> rootedVerticesOld{ExternalPolyakovLoopOld{1,true}, ExternalPolyakovLoopOld{2,false}};

    PureGaugeWeight<double> someTest(WeightContainer, rootedVertices);
    PureGaugeWeight<GiNaC::numeric> anotherTest(WeightContainer, rootedVertices);
    PureGaugeWeightOld yetAnotherTest(WeightContainer, rootedVerticesOld);
    //StaticQuarkWeight staticQuarkTest(WeightContainer, WeightContainer, rootedVertices);

    //std::cout << staticQuarkTest.Weight() << "\n";

    g6String = "A_"; /// one-link graph
    tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer OneLinkContainer(2, m, g); /// container from densenauty
    /// vertices 1 and 2 are rooted
    //WeightContainer.SetRootedVertex(0,0);
    //WeightContainer.SetRootedVertex(1,1);
    std::cout << OneLinkContainer;

    StaticQuarkWeight staticQuarkOneLink(OneLinkContainer, OneLinkContainer, std::vector<ExternalPolyakovLoop>());

    auto tempWeight = staticQuarkOneLink.Weight();
    //std::cout << tempWeight << "\n";

    return 0;

    std::cout << someTest.Weight() << " " << anotherTest.Weight() << " " << yetAnotherTest.Weight() << "\n";

    CubicLattice MyLattice(100);
    std::vector<unsigned int> indicesCenter(MyLattice.GetDim(), MyLattice.GetN()/2); /// spatial indices for first vertex
    auto index1 = MyLattice.GetSiteIndex(indicesCenter); /// index for vertex 1
    auto index2 = MyLattice.GetNearestNeighbor(index1, 1);
    auto index3 = MyLattice.GetNearestNeighbor(index1, 2);
    auto index4 = MyLattice.GetNearestNeighbor(index2, 2);

    VertexEmbedList myEmbedList(MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor);
    myEmbedList.AddFixedVertexEmbed(0 /*fixedNbr (color)*/, 1 /* vertexNbr (always 1 after NAUTY does canonicalization) */, index1);
    myEmbedList.AddFixedVertexEmbed(1 /*fixedNbr (color)*/, 2 /* vertexNbr (always 2 after NAUTY does canonicalization) */, index2);
    myEmbedList.AddVertexEmbed(3, index3);
    myEmbedList.AddVertexEmbed(4, index4);

    ZClusterPureGaugeArbEmbedding myZ(WeightContainer, myEmbedList, &MyLattice, std::vector<bool>{true,false});

    std::cout << "TESTZCLUSTER: " << myZ.ComputeLambdaPolynomial().GetPolynomial() << "\n";

    TwoPointCorrelator<ZClusterPureGaugeArbEmbedding> myCorr(WeightContainer, myEmbedList, &MyLattice);
    std::cout << "TESTFULLCORR:" << myCorr.GetFullCorrelatorGiNaC() << "\n";
    std::cout << "TESTEXPANDEDCORR:" << myCorr.GetExpandedCorrelatorGiNaC() << "\n";

    CanonicalGraphManager MyManager(5);

    XiRecursionRooted<ZClusterPureGaugeArbEmbedding> myRootedRecursion(&MyManager, WeightContainer, myEmbedList, &MyLattice);

    auto fullXi = myRootedRecursion.GetFullXiGiNaC();
    std::cout << "FULL_XI: " << fullXi << "\n";

    auto expandedXi = myRootedRecursion.GetExpandedXiGiNaC();
    std::cout << "EXPANDED_XI: " << expandedXi << "\n";

    return 0;
}
