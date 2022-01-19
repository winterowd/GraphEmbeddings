#include "ZClusterPureGaugeArbEmbedding.h"

ZClusterPureGaugeArbEmbedding::ZClusterPureGaugeArbEmbedding(GraphContainer *container, VertexEmbedList *clusterEmbedList, CubicLattice* lattice) :
    ClusterContainer(container),
    ClusterEmbedList(clusterEmbedList),
    Lattice(lattice)
{
    /// set up variables by identifying each bond of embedding
    auto tempOneLink = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreNN(index1, index2); };
    this->FillBondType(this->OneLink, tempOneLink, "ONE_LINK");
    this->AllTotalBondCountsPlusOne[0] = this->OneLink.size()+1; /// set size
    if (this->OneLink.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[1];

    auto tempSquareDiagonal = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreNNN(index1, index2); };
    this->FillBondType(this->SquareDiagonal, tempSquareDiagonal, "SQUARE_DIAGONAL");
    this->AllTotalBondCountsPlusOne[1] = this->SquareDiagonal.size()+1; /// set size
    if (this->SquareDiagonal.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[1];

    auto tempStraightTwoLink = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreThirdNN(index1, index2); };
    this->FillBondType(this->StraightTwoLink, tempStraightTwoLink, "STRAIGHT_TWO_LINK");
    this->AllTotalBondCountsPlusOne[2] = this->StraightTwoLink.size()+1; /// set size
    if (this->StraightTwoLink.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[2];

    auto tempCubeDiagonal = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreFourthNN(index1, index2); };
    this->FillBondType(this->CubeDiagonal, tempCubeDiagonal, "CUBE_DIAGONAL");
    this->AllTotalBondCountsPlusOne[3] = this->CubeDiagonal.size()+1; /// set size
    if (this->CubeDiagonal.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[3];

    this->TotalBondCounts = this->OneLink.size()+this->SquareDiagonal.size()+this->StraightTwoLink.size()+this->CubeDiagonal.size(); /// precompute total number of bonds (length of valid permutation)
    this->ZCoefficients.resize(this->LinearIndexMax); /// resize
    this->EvaluateZ(); /// evaluate all terms
}

/// convert linear index to a tuple corresponding to the powers of each coupling \lambda_i
/// @param index: linear index
std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings> ZClusterPureGaugeArbEmbedding::LinearIndexToPowersOfCouplings(int index) const
{
#ifdef DEBUG
    if (index < 0 || index >= this->LinearIndexMax)
        throw std::invalid_argument("ERROR: LinearIndexToPowersOfCouplings requires 0<= index < NTotalBondCounts!\n");
#endif
    std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings> result{0, 0, 0, 0};
    int temp = index;
    for (int i=0; i<this->NbrCouplings; ++i)
    {
        if (this->AllTotalBondCountsPlusOne[i]!=0)
        {
            result[i] = temp%this->AllTotalBondCountsPlusOne[i];
            temp = temp/this->AllTotalBondCountsPlusOne[i];
        }
    }
    return result;
}

/// check if a given tuple corresponding to powers of couplings is valid
/// @param powers: array containing powers of the couplings
bool ZClusterPureGaugeArbEmbedding::ValidPowersOfCouplings(const std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>& powers) const
{
    for (int i=0; i<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++i)
        if (powers[i] < 0 || powers[i] >= this->AllTotalBondCountsPlusOne[i])
            return false;
    return true;
}

/// convert from a tuple corresponding to powers of couplings to a linear index
/// @param powers: powers of the couplings \lambda_i
int ZClusterPureGaugeArbEmbedding::PowersOfCouplingsToLinearIndex(const std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>& powers) const
{
    if (!this->ValidPowersOfCouplings(powers))
        throw std::invalid_argument("ERROR: PowersOfCouplingsToLinearIndex requires powers to contain integers in the appropriate range!\n");

    int linearIndex = 0;
    int step = 1;
    for (int i=0; i<this->NbrCouplings; ++i)
    {
        linearIndex += powers[i]*step;
        step *= this->AllTotalBondCountsPlusOne[i];
    }
    return linearIndex;
}

/// fill the vector with certain type of edges (defined by isEdgeValid)
void ZClusterPureGaugeArbEmbedding::FillBondType(std::vector<UndirectedEdge>& result, std::function<bool(unsigned int, unsigned int)> isEdgeValid, std::string bondType)
{
    for (int i=0; i<this->ClusterContainer->GetL(); ++i)
    {
        auto tempEdge = this->ClusterContainer->GetEdge(i);
        if (isEdgeValid(this->ClusterEmbedList->GetVertexSiteIndex(tempEdge.FirstVertex), this->ClusterEmbedList->GetVertexSiteIndex(tempEdge.SecondVertex)))
        {
#ifdef DEBUG
            std::vector<unsigned int> indices1(3), indices2(3);
            this->Lattice->GetSiteCoordinates(this->ClusterEmbedList->GetVertexSiteIndex(tempEdge.FirstVertex), indices1);
            this->Lattice->GetSiteCoordinates(this->ClusterEmbedList->GetVertexSiteIndex(tempEdge.SecondVertex), indices2);
            std::cout << "bond number " << i+1 << " found to be of type " << bondType << "!\n";
            std::cout << "Between vertex " << tempEdge.FirstVertex << " at ( ";
            for (int j=0; j<3; ++j)
                std::cout << indices1[j] << ", ";
            std::cout << ") and vertex " << tempEdge.SecondVertex << " at ( ";
            for (int j=0; j<3; ++j)
                std::cout << indices2[j] << ", ";
            std::cout << ")!\n";
#endif
            result.push_back(tempEdge);
        }
    }
#ifdef DEBUG
    std::cout << "Identified " << result.size() << " bonds of type " << bondType << "!\n";
#endif
}

/// recursively generate all combinations of boolean string
/// @param tmp: current string
/// @param nbrBondsRemaining: number of bonds/booleans to add to current string
void ZClusterPureGaugeArbEmbedding::GenerateTermsIntegrand(std::vector<bool>& tmp, int nbrBondsRemaining)
{
    if (nbrBondsRemaining==0)
    {
#ifdef DEBUG
        if (tmp.size()!=this->TotalBondCounts)
            throw std::invalid_argument("ERROR: GenerateTermsIntegrand only add permutations to IntegrandTerms if they are of size TotalBondCounts!\n");
#endif
        this->IntegrandTerms.push_back(tmp);
        return;
    }

    /// true: bond exists ( \lambda_i (L_x L^*_y + L^*_x L_y), i=1,2,3 )
    tmp.push_back(true);
    this->GenerateTermsIntegrand(tmp, nbrBondsRemaining-1);
    tmp.pop_back();

    /// false: bond does not exist (1)
    tmp.push_back(false);
    this->GenerateTermsIntegrand(tmp, nbrBondsRemaining-1);
    tmp.pop_back();

    return;
}

/// evaluate all of the terms in the partition function
/// store results in ZCoefficients
void ZClusterPureGaugeArbEmbedding::EvaluateZ()
{
    std::vector<bool> tmp;
    this->GenerateTermsIntegrand(tmp, this->TotalBondCounts);
    for (int i=0; i<this->IntegrandTerms.size(); ++i)
    {
        auto edgesAndBondCount = this->ZIntegrandToUndirectedEdgesAndBondCounts(this->IntegrandTerms[i]); /// get edges and bond counts
        auto relabeledEdgesAndVertexMap = SubDiagramGenerator::GetRelabeledEdgesAndVertexMap(edgesAndBondCount.first); /// relabeled edges and vertex map
        GraphContainer tempContainer(relabeledEdgesAndVertexMap.second.size(), 1, relabeledEdgesAndVertexMap.first); /// graph container
        PureGaugeWeight tempWeightObject(&tempContainer);
        double weight = tempWeightObject.Weight(); /// calculate weight of graph
#ifdef DEBUG
        //// print out i, bond counts, Container, VertexMap
        std::cout << "**** DEBUG_EVALUATEZ ****\n";
        std::cout << "term " << i << " with bond_count:";
        for (int j=0; j<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++j)
            std::cout << " " << edgesAndBondCount.second[j];
        std::cout << "\n";
        std::cout << tempContainer;
        for (int j=0; j<relabeledEdgesAndVertexMap.second.size(); ++j)
            std::cout << "Vertex " << j+1 << " maps to original vertex " << relabeledEdgesAndVertexMap.second[j] << "\n";
        std::cout << "Weight: " << weight << "\n";
        std::cout << "**** DEBUG_EVALUATEZ ****\n";
#endif
        this->ZCoefficients[this->PowersOfCouplingsToLinearIndex(edgesAndBondCount.second)] += weight;
    }
}

/// convert the boolean string to undirected edges and counts for each bond type
/// @param permutation: boolean string of length \sum_i N_i
std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>> ZClusterPureGaugeArbEmbedding::ZIntegrandToUndirectedEdgesAndBondCounts(const std::vector<bool>& permutation)
{
#ifdef DEBUG
    if (this->TotalBondCounts!=permutation.size())
        throw std::invalid_argument("ERROR: ZIntegrandToUndirectedEdges requires permutation to be of size TotalBondCounts!\n");
#endif
    std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings> bondCounts{};
    std::vector<UndirectedEdge> result;
    for (int i=0; i<this->OneLink.size(); ++i) /// translate one-link to edges
    {
        if (permutation[i]) /// is this bond requested?
        {
            bondCounts[0]++;
            result.push_back(this->OneLink[i]);
        }
    }
    int offset = this->OneLink.size(); /// offset for accessing elements of permutation
    for (int i=0; i<this->SquareDiagonal.size(); ++i) /// translate hook to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[1]++;
            result.push_back(this->SquareDiagonal[i]);
        }
    }
    offset += this->SquareDiagonal.size(); /// increment offset
    for (int i=0; i<this->StraightTwoLink.size(); ++i) /// translate straight two-link to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[2]++;
            result.push_back(this->StraightTwoLink[i]);
        }
    }
    offset += this->StraightTwoLink.size(); /// increment offset
    for (int i=0; i<this->CubeDiagonal.size(); ++i) /// translate cube diagonals to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[3]++;
            result.push_back(this->CubeDiagonal[i]);
        }
    }
#ifdef DEBUG
    for (int i=0; i<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++i)
        if (bondCounts[i]>this->AllTotalBondCountsPlusOne[i])
            throw std::invalid_argument("ERROR: ZIntegrandToUndirectedEdges requires bondCounts to be less than corresponding element of AllTotalBondCountsPlusOne!\n");
#endif
    return std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>>(result, bondCounts);
}

/// print out all non-zero coefficients of Z
void ZClusterPureGaugeArbEmbedding::PrintZ() const
{
    std::cout << "**** PrintZ ****\n";
    for (int i=0; i<this->LinearIndexMax; ++i)
    {
        double tempCoefficient = this->ZCoefficients[i];
        if (tempCoefficient!=0)
        {
            auto powersCouplings = this->LinearIndexToPowersOfCouplings(i);
            std::cout << "( ";
            for (int j=0; j<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++j)
                std::cout << powersCouplings[j] << ", ";
            std::cout << "): " << tempCoefficient << "\n";
        }
    }
    std::cout << "**** PrintZ ****\n";
}

/// print the contributions at a given order in the couplings
/// @param powers: requested powers of each coupling \lambda_i
void ZClusterPureGaugeArbEmbedding::PrintContributionZFixedOrder(const std::array<int, ZClusterPureGaugeArbEmbedding::NbrCouplings>& powers)
{
    if (!this->ValidPowersOfCouplings(powers))
        throw std::invalid_argument("ERROR: PrintContributionZFixedOrder requires powers to contain integers in the appropriate range!\n");

    std::cout << "**** PrintContributionZFixedOrder ****\n";
    std::cout << "( ";
    for (int j=0; j<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++j)
        std::cout << powers[j] << ", ";
    std::cout << ")\n";
    for (int i=0; i<this->IntegrandTerms.size(); ++i)
    {
        auto edgesAndBondCount = this->ZIntegrandToUndirectedEdgesAndBondCounts(this->IntegrandTerms[i]); /// get edges and bond counts
        bool correctBondCount = true;
        for (int j=0; j<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++j)
            if (powers[j]!=edgesAndBondCount.second[j])
                correctBondCount = false;
        if (correctBondCount)
        {
            auto relabeledEdgesAndVertexMap = SubDiagramGenerator::GetRelabeledEdgesAndVertexMap(edgesAndBondCount.first); /// relabeled edges and vertex map
            GraphContainer tempContainer(relabeledEdgesAndVertexMap.second.size(), 1, relabeledEdgesAndVertexMap.first); /// graph container
            PureGaugeWeight tempWeightObject(&tempContainer);
            double weight = tempWeightObject.Weight(); /// calculate weight of graph
            std::cout << "term " << i << " with bond_count:";
            for (int j=0; j<ZClusterPureGaugeArbEmbedding::NbrCouplings; ++j)
                std::cout << " " << edgesAndBondCount.second[j];
            std::cout << "\n";
            std::cout << tempContainer;
            for (int j=0; j<relabeledEdgesAndVertexMap.second.size(); ++j)
                std::cout << "Vertex " << j+1 << " maps to original vertex " << relabeledEdgesAndVertexMap.second[j] << "\n";
            std::cout << "Weight: " << weight << "\n";
        }
    }
    std::cout << "**** PrintContributionZFixedOrder ****\n";
}
