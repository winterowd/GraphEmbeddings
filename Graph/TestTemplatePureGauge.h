#ifndef TESTTEMPLATEPUREGAUGE_H
#define TESTTEMPLATEPUREGAUGE_H

#include <ginac/ginac.h>

#include "GraphContainer.h"
#include "PureGaugeweight.h"

template <typename T>
class TestTemplatePureGauge
{
private:
    GraphContainer Container; /// pointer to graph container

    std::vector<UndirectedEdge> Edges; /// list of (undirected) edges

    std::vector<ExternalPolyakovLoop> ExternalVertices; /// list of vertices where external L or L* is placed

    T TotalWeight; /// total weight of graph using all possible combinations of directed bonds

    void GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining);

    T GetGraphWeightFixedBonds(const std::vector<bool> &directedEdges);

    std::vector<SiteCount> GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges);

    double BinomialCoefficient(int n, int k);

    double Factorial(int n);

    T SingleSiteIntegral(int n, int m);


public:
    TestTemplatePureGauge(const GraphContainer& container, const std::vector<ExternalPolyakovLoop>& externalVertices);

    T Weight(); /// routine called by user (interface)
};

#endif // TESTTEMPLATEPUREGAUGE_H
