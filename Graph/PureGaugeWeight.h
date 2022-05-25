#ifndef PUREGAUGEWEIGHT_H
#define PUREGAUGEWEIGHT_H

#include <ginac/ginac.h>

#include "GraphContainer.h"
#include "ExternalPolyakovLoop.h"

template <typename T>
class PureGaugeWeight
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
    PureGaugeWeight(const GraphContainer& container, const std::vector<ExternalPolyakovLoop>& externalVertices);

    T Weight(); /// routine called by user (interface)
};

#endif // PUREGAUGEWEIGHT_H
