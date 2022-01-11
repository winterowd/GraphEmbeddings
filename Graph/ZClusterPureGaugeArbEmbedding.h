#ifndef ZCLUSTERPUREGAUGEARBEMBEDDING_H
#define ZCLUSTERPUREGAUGEARBEMBEDDING_H

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "CubicLattice.h"
#include "SubDiagramGenerator.h"
#include "PureGaugeweight.h"

class ZClusterPureGaugeArbEmbedding
{
private:
    /**** private variables ****/
    GraphContainer *ClusterContainer; /// container corresponding to the cluster

    VertexEmbedList *ClusterEmbedList; /// embed list corresponding to the cluster

    CubicLattice *Lattice; /// lattice object

    /// all permutations of a string of bools of length \sum_i N_i
    std::vector<std::vector<bool>> IntegrandTerms; /// terms in the integrand of Z expressed as strings of booleans

    /// coefficient of each term  \prod_i \lambda^{n_i}_i (linear index)
    std::vector<double> ZCoefficients; /// size \prod_i (N_i+1) (where N_i are the number of each type of bond for cluster)

    /// lists of bonds
    std::vector<UndirectedEdge> OneLink; /// one-link edges

    std::vector<UndirectedEdge> SquareDiagonal; /// "square"-diagonal edges (of length \sqrt{2})

    std::vector<UndirectedEdge> StraightTwoLink; /// straight two-link edges

    std::vector<UndirectedEdge> CubeDiagonal; /// "cube"-diagonal edges (of length \sqrt{3})

public:
    static const int NbrCouplings = 4; /// NN, NNN, Third-NN, Fourth-NN
private:
    std::array<int, NbrCouplings> AllTotalBondCountsPlusOne; /// N_i+1, i=1,..,NbrCouplings

    int TotalBondCounts; /// sum of number of OneLink, Hooks, and StraightTwoLink

    int LinearIndexMax; /// \prod_i (N_i+1), where N_i is the total available bonds of a given type

    /**** private methods ****/
    std::array<int, NbrCouplings> LinearIndexToPowersOfCouplings(int index) const; /// convert from linear index to powers of couplings

    int PowersOfCouplingsToLinearIndex(const std::array<int, NbrCouplings>& powers) const; /// convert from powers of couplings to linear index

    bool ValidPowersOfCouplings(const std::array<int, NbrCouplings>& powers) const; /// check if couplings are in valid range

    /// get all edges which correspond to two-point interactions of a given type, \lambda_i
    void FillBondType(std::vector<UndirectedEdge>& result, std::function<bool(unsigned int, unsigned int)> isEdgeValid, std::string bondType);

    /// convert a permutation string to edges and powers of the couplings i.e. bond counts
    std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>> ZIntegrandToUndirectedEdgesAndBondCounts(const std::vector<bool>& permutation);

    /// generate all "terms" i.e. permutations recursively
    void GenerateTermsIntegrand(std::vector<bool>& tmp, int nbrBondsRemaining);

    /// evalute the partition function
    void EvaluateZ();

public:
    ZClusterPureGaugeArbEmbedding(GraphContainer *container, VertexEmbedList *clusterEmbedList, CubicLattice* lattice);

    /**** accessors ****/
    int GetNbrSquareDiagonal() const { return this->SquareDiagonal.size(); }

    int GetNbrOneLink() const { return this->OneLink.size(); }

    int GetNbrStraightTwoLink() const { return this->StraightTwoLink.size(); }

    int GetNbrCubeDiagonal() const { return this->CubeDiagonal.size(); }

    /**** debugging routines ****/
    void PrintZ() const;

    void PrintContributionZFixedOrder(const std::array<int, NbrCouplings>& powers);
};

#endif // ZCLUSTERPUREGAUGEARBEMBEDDING_H
