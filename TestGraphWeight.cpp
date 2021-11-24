#include <iostream>
#include <chrono>

#include "nauty.h"
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

    PureGaugeWeight MyWeightObject(&WeightContainer);

    auto result = MyWeightObject.Weight();
    std::cout << "TOTAL_WEIGHT: " << result << "\n";

    DYNFREE(g, g_sz); /// free graph

    n = 5;
    m = SETWORDSNEEDED(n);

    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string inputFilename = "graphs_g6_connected_order_5.dat";

    FILE *fp = fopen(inputFilename.c_str(), "r"); /// open file containing graphs of a fixed order from which we will generate two-rooted graphs
    if (fp!=NULL)
        std::cout << "TestGraphWeight: Successfully opened " << inputFilename << std::endl;
    else
        throw std::invalid_argument("TestGraphWeight: Error opening "+inputFilename);

    GraphContainer refContainer(n, m);
    while (1) /// read graphs in file (ASSUME THEY ARE ALL OF THE SAME ORDER N!)
    {
        graph *g1 = readg(fp,g,0,&m,&n);
        if (g1 == NULL)
            break;

        refContainer.SetGraphFromDenseNauty(g);
        refContainer.PrintM();
        PureGaugeWeight weight(&refContainer);
        result = weight.Weight();
        std::cout << "TOTAL_WEIGHT: " << result << "\n";

    }

    DYNFREE(g, g_sz); /// free graph

}
