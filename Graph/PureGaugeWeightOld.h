#ifndef PUREGAUGEWEIGHTOLD_H
#define PUREGAUGEWEIGHTOLD_H

#include "GraphContainer.h"
#include <math.h>
#include <iostream>

struct ExternalPolyakovLoopOld {
    int Label; /// vertex label
    bool Fundamental; /// fundamental (L) or antifundamental (L*)?
};

class PureGaugeWeightOld
{
private:
    GraphContainer Container; /// pointer to graph container

    std::vector<UndirectedEdge> Edges; /// list of (undirected) edges

    std::vector<ExternalPolyakovLoopOld> ExternalVertices; /// list of vertices where external L or L* is placed

    double TotalWeight; /// total weight of graph using all possible combinations of directed bonds

    double GetGraphWeightFixedBonds(const std::vector<bool> &directedEdges);

    std::vector<SiteCount> GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges);

    double SingleSiteIntegral(int n, int m); /// single site integration

    void GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining);

    double BinomialCoefficient(int n, int k);

    double Factorial(int n);

public:
    PureGaugeWeightOld(const GraphContainer& container, const std::vector<ExternalPolyakovLoopOld>& externalVertices);

    double Weight(); /// routine called by user (interface)
};

#endif // PUREGAUGEWEIGHTOLD_H
