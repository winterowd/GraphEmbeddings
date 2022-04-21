#ifndef TWOPOINTCORRELATOR_H
#define TWOPOINTCORRELATOR_H

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "ZClusterPureGaugeArbEmbedding.h"
#include "CubicLattice.h"

class TwoPointCorrelator
{
private:
    GraphContainer ContainerRootedCluster; /// container for two-rooted graph

    VertexEmbedList EmbedListRootedCluster; /// embed list for two-rooted graph

    CubicLattice *Lattice; /// pointer to lattice object

    int MaxManhattanDistance; /// maximum Manhattan distance

    /// pointers to cluster objects (dynamically allocate as they do not have default constructor)
    std::vector<ZClusterPureGaugeArbEmbedding<GiNaC::numeric>> CorrTerms;

public:
    TwoPointCorrelator(const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice* lattice, int maxManhattanDistance=10);

    GiNaC::ex GetFullCorrelatorGiNaC();

    GiNaC::ex GetExpandedCorrelatorGiNaC();

    void PrintCorrelatorTerms(); /// debugging routine
};

#endif // TWOPOINTCORRELATOR_H
