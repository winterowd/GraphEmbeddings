#ifndef STATICQUARKWEIGHT_H
#define STATICQUARKWEIGHT_H

#include "GraphContainer.h"
#include "ExternalPolyakovLoop.h"
#include "AuxiliaryRoutinesForGINAC.h"

#include <ginac/ginac.h>

/// Single flavor with static quarks
class StaticQuarkWeight
{
private:
    GraphContainer Container; /// pointer to graph container (g)

    std::vector<UndirectedEdge> Edges; /// list of (undirected) edges

    std::vector<ExternalPolyakovLoop> ExternalVertices; /// list of vertices where external L or L* is placed

    GiNaC::ex TotalWeight; /// total weight of graph using all possible combinations of directed bonds

    static bool StaticDeterminantPrecomputed; /// flag

    void GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining);

    GiNaC::ex GetGraphWeightFixedBonds(const std::vector<bool> &directedEdges);

    std::vector<SiteCount> GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges);

    GiNaC::ex ComputeSiteContribution(int nbrIn, int nbrOut);

    GiNaC::numeric SingleSiteIntegral(int n, int m);

public:
    StaticQuarkWeight(const GraphContainer& container, const std::vector<ExternalPolyakovLoop>& externalVertices);

    GiNaC::ex Weight(); /// routine called by user (interface)

    static std::vector<std::pair<GiNaC::ex, SiteCount>> StaticDeterminant;

    static GiNaC::ex SingleSiteIntegratedStaticDeterminant;

};

// some useful functions
std::vector<std::pair<GiNaC::ex, SiteCount>> ComputeStaticDeterminant();
GiNaC::ex ComputeSingleSiteIntegratedStaticDeterminant();

#endif // STATICQUARKWEIGHT_H
