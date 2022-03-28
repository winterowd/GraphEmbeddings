#include "AuxiliaryRoutinesForGINAC.h"
#include "MyLambdaPolynomial.h"
#include "PureGaugeweight.h"
#include "PureGaugeWeightOld.h"

int main()
{
    std::vector<std::pair<double, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients1;
    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients2;
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder1{1,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder2{3,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder3{1,1,0,0};

    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1., tempOrder1});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{3., tempOrder2});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1./6., tempOrder3});

    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder1});
    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(3), tempOrder2});
    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1,6), tempOrder3});

    auto result1 = AuxiliaryRoutinesForGinac::PolynomialFromCoefficients<double>(coefficients1, 6);
    auto result2 = AuxiliaryRoutinesForGinac::PolynomialFromCoefficients<GiNaC::numeric>(coefficients2, 6);

    std::cout << "RESULT1: " << result1 << "\n";
    std::cout << "RESULT2: " << result2 << "\n";

    auto result3 = MyLambdaPolynomial<double>(coefficients1, 6);

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

    GraphContainer WeightContainer(n, m, g); /// container from densenauty
    std::vector<ExternalPolyakovLoop> rootedVertices{ExternalPolyakovLoop{1,true}, ExternalPolyakovLoop{2,false}};
    std::vector<ExternalPolyakovLoopOld> rootedVerticesOld{ExternalPolyakovLoopOld{1,true}, ExternalPolyakovLoopOld{2,false}};

    PureGaugeWeight<double> someTest(WeightContainer, rootedVertices);
    PureGaugeWeight<GiNaC::numeric> anotherTest(WeightContainer, rootedVertices);
    PureGaugeWeightOld yetAnotherTest(WeightContainer, rootedVerticesOld);

    std::cout << someTest.Weight() << " " << anotherTest.Weight() << " " << "\n"; //yetAnotherTest.Weight() << "\n";

    return 0;
}
