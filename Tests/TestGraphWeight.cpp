#include <iostream>
#include <chrono>

#include "nauty.h"

#include "PureGaugeWeightOld.h"
#include "PureGaugeweight.h"

int main(int argc, char *argv[])
{
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

    GraphContainer WeightContainer(n, m, g); /// container from densenauty

    std::vector<ExternalPolyakovLoopOld> rootedVerticesOld1{ExternalPolyakovLoopOld{1,true}, ExternalPolyakovLoopOld{2,false}};
    PureGaugeWeightOld MyWeightObjectOld1(WeightContainer, rootedVerticesOld1);
    auto resultOld = MyWeightObjectOld1.Weight();
    std::cout << "TOTAL_WEIGHTOLD1: " << resultOld << "\n";

    std::vector<ExternalPolyakovLoop> rootedVertices1{ExternalPolyakovLoop{1,true}, ExternalPolyakovLoop{2,false}};
    PureGaugeWeight<GiNaC::numeric> MyWeightObject1(WeightContainer, rootedVertices1);
    auto result = MyWeightObject1.Weight();
    std::cout << "TOTAL_WEIGHT1: " << result << "\n";

    std::vector<ExternalPolyakovLoopOld> rootedVerticesOld2{ExternalPolyakovLoopOld{1,true}};
    PureGaugeWeightOld MyWeightObjectOld2(WeightContainer, rootedVerticesOld2);
    resultOld = MyWeightObjectOld2.Weight();
    std::cout << "TOTAL_WEIGHTOLD2: " << resultOld << "\n";

    std::vector<ExternalPolyakovLoop> rootedVertices2{ExternalPolyakovLoop{1,true}};
    PureGaugeWeight<GiNaC::numeric> MyWeightObject2(WeightContainer, rootedVertices2);
    result = MyWeightObject2.Weight();
    std::cout << "TOTAL_WEIGHT2: " << result << "\n";

    std::vector<ExternalPolyakovLoopOld> rootedVerticesOld3{ExternalPolyakovLoopOld{2,false}};
    PureGaugeWeightOld MyWeightObjectOld3(WeightContainer, rootedVerticesOld3);
    resultOld = MyWeightObjectOld3.Weight();
    std::cout << "TOTAL_WEIGHTOLD3: " << resultOld << "\n";

    std::vector<ExternalPolyakovLoop> rootedVertices3{ExternalPolyakovLoop{2,false}};
    PureGaugeWeight<GiNaC::numeric> MyWeightObject3(WeightContainer, rootedVertices3);
    result = MyWeightObject3.Weight();
    std::cout << "TOTAL_WEIGHT3: " << result << "\n";

    DYNFREE(g, g_sz); /// free graph

}
