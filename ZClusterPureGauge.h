#ifndef ZCLUSTERPUREGAUGE_H
#define ZCLUSTERPUREGAUGE_H

#include <vector>
#include <iostream>

#include "GraphContainer.h"
#include "VertexEmbedList.h"
#include "SubDiagramGenerator.h"
#include "PureGaugeweight.h"

/// Helper class which stores identified two-point interactions from subdiagrams of a given graph G
/// used ONLY in conjunction with SubDiagramGenerator object and a fixed graph embedding of a given graph G
class InteractionSubDiagram
{
private:
    std::vector<int> Ends; /// vertex numbers

    int SubDiagramLabel; /// label for corresponding subdiagram (access through SubDiagramGenerator object)
public:
    /// constructor
    InteractionSubDiagram(int v1, int v2, int label) :
        Ends(2),
        SubDiagramLabel(label)
    {
        if (v1==v2)
            throw std::invalid_argument("ERROR: InteractionSubDiagram requires unique vertices v1 and v2!\n");
        if (v1<v2) /// order for comparison
        {
            this->Ends[0] = v1;
            this->Ends[1] = v2;
        }
        else
        {
            this->Ends[0] = v2;
            this->Ends[1] = v1;
        }
    }

    /// accessors
    int GetEnd(int index) const { if (index!=0 && index!=1) throw std::invalid_argument("ERROR: GetEnd expects index to be 0 or 1!\n"); return this->Ends[index]; }
    int GetSubDiagramLabel() const { return this->SubDiagramLabel; }

    /// relational opeartors
    friend bool operator==(const InteractionSubDiagram& lhs, const InteractionSubDiagram& rhs);
    friend bool operator<(const InteractionSubDiagram& lhs, const InteractionSubDiagram& rhs);

};

/// equality based on Ends only (sorted)
inline bool operator==(const InteractionSubDiagram& lhs, const InteractionSubDiagram& rhs)
{
    return (lhs.Ends[0]==rhs.Ends[0] && lhs.Ends[1]==rhs.Ends[1]);
}

/// inequality operator
inline bool operator!=(const InteractionSubDiagram& lhs, const InteractionSubDiagram& rhs)
{
    return !(rhs==lhs);
}

/// lexicographical order!
inline bool operator<(const InteractionSubDiagram& lhs, const InteractionSubDiagram& rhs)
{
    for (int i=0; i<2; ++i)
    {
        if (lhs.Ends[i]!=rhs.Ends[i])
            return lhs.Ends[i]<rhs.Ends[i];
    }
    return false;
}

class ZClusterPureGauge
{
private:
    /**** private variables ****/
    GraphContainer *ClusterContainer; /// container corresponding to the cluster

    VertexEmbedList *ClusterEmbedList; /// embed list corresponding to the cluster

    CubicLattice *Lattice; /// lattice object

    SubDiagramGenerator MySubDiagramGenerator; /// subdiagrams of cluster

    /// all permutations of a string of bools of length \sum_i N_i
    std::vector<std::vector<bool>> IntegrandTerms; /// terms in the integrand of Z expressed as strings of booleans

    /// coefficient of each term  \prod_i \lambda^{n_i}_i (linear index)
    std::vector<double> ZCoefficients; /// size \prod_i (N_i+1) (where N_i are the number of each type of bond for cluster)

    /// lists of relevant subgraphs
    std::vector<InteractionSubDiagram> OneLink; /// indices corresponding to one-link subgraphs

    std::vector<InteractionSubDiagram> Hooks; /// indices corresponding to UNIQUE hooks

    std::vector<InteractionSubDiagram> StraightTwoLink; /// indices corresponding to straight two-link subgraphs

    std::vector<InteractionSubDiagram> CubeDiagonal; /// indices corresponding to UNIQUE diagonals of a cube

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

    void FillOneLink(); /// all one-link subdiagrams

    /// get all subdiagrams whose ENDPOINTS correspond to two-point interactions of a given type, \lambda_i
    void FillBondType(std::vector<InteractionSubDiagram>& result, std::function<bool(unsigned int, unsigned int)> isValidNeighbor, int validNbrVertices, int validNbrBonds, std::string bondType);

    /// check pairs of verices in a VertexEmbedList and see if they are separated by a distance corresponding to \lambda_i (assume there is at most ONE PAIR of vertices satisfying this)
    bool IsSubDiagramValidNhbrPath(const VertexEmbedList& list, std::vector<int>& nnnLabels, std::function<bool(unsigned int, unsigned int)> isValidNeighbor, int validNbrVertices);

    /// convert a permutation string to edges and powers of the couplings i.e. bond counts
    std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGauge::NbrCouplings>> ZIntegrandToUndirectedEdgesAndBondCounts(const std::vector<bool>& permutation);

    /// generate all "terms" i.e. permutations recursively
    void GenerateTermsIntegrand(std::vector<bool>& tmp, int nbrBondsRemaining);

    /// evalute the partition function
    void EvaluateZ();

public:
    ZClusterPureGauge(GraphContainer* clusterContainer, VertexEmbedList* clusterEmbedList, CubicLattice* lattice); /// constructor

    /**** accessors ****/
    int GetNbrHooks() const { return this->Hooks.size(); }

    int GetNbrOneLink() const { return this->OneLink.size(); }

    int GetNbrStraightTwoLink() const { return this->StraightTwoLink.size(); }

    /**** debugging routines ****/
    void PrintZ() const;

    void PrintContributionZFixedOrder(const std::array<int, NbrCouplings>& powers);

};

#endif // ZCLUSTERPUREGAUGE_H
