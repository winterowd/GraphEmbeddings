#ifndef TWOPOINTCORRELATOR_H
#define TWOPOINTCORRELATOR_H

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "ZClusterPureGaugeArbEmbedding.h"
#include "ZClusterStaticQuarkArbEmbedding.h"
#include "CubicLattice.h"

template<typename T>
class TwoPointCorrelator
{
private:
    GraphContainer ContainerRootedCluster; /// container for two-rooted graph

    VertexEmbedList EmbedListRootedCluster; /// embed list for two-rooted graph

    CubicLattice *Lattice; /// pointer to lattice object

    int MaxManhattanDistance; /// maximum Manhattan distance

    int MaxOrderH1; /// maximum order in h_1

    int MaxOrderHBar1; /// maximum order in \bar{h}_1

    /// vector of cluster objects (dynamically allocate as they do not have default constructor)
    std::vector<T> CorrTerms;

public:
    TwoPointCorrelator(const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice* lattice, int maxManhattanDistance=10, int maxOrderH1=0, int maxOrderHBar1=0);

    GiNaC::ex GetFullCorrelatorGiNaC();

    GiNaC::ex GetExpandedCorrelatorGiNaC();

    int GetMaxOrderH1() const { return this->MaxOrderH1; }

    int GetMaxOrderHBar1() const { return this->MaxOrderHBar1; }

    void PrintCorrelatorTerms(); /// debugging routine
};

#endif // TWOPOINTCORRELATOR_H
