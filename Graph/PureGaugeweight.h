#ifndef PUREGAUGEWEIGHT_H
#define PUREGAUGEWEIGHT_H

#include "GraphContainer.h"
#include <math.h>
#include <iostream>

class PureGaugeWeight
{
private:
    GraphContainer *Container; /// pointer to graph container

    std::vector<UndirectedEdge> Edges; /// list of (undirected) edges

    double TotalWeight; /// total weight of graph using all possible combinations of directed bonds

    double GetGraphWeightFixedBonds(const std::vector<bool> &directedEdges);

    std::vector<SiteCount> GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges);

    double SingleSiteIntegral(int n, int m); /// single site integration

    void GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining);

    double BinomialCoefficient(int n, int k);

    double Factorial(int n);

public:
    PureGaugeWeight(GraphContainer *container);

    double Weight(); /// routine called by user (interface)
};

#endif // PUREGAUGEWEIGHT_H
