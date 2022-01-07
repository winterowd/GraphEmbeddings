#include "ZClusterPureGauge.h"

/// constructor: takes a graph, a VALID embedding on the cubic lattice, and a lattice object.
/// Creates lists of all subdiagrams of cluster which can be mapped to a unique bond of ith order (\lambda_i)
/// Evaluates Z and saves the coefficients
/// @param clusterContainer: pointer to a GraphContainer object for G (memory managed externally)
/// @param clusterEmbedList: pointer to a VertexEmbedList object for G on the cubic lattice (memory managed externally)
/// @param lattice: pointer to CubicLattice object
ZClusterPureGauge::ZClusterPureGauge(GraphContainer *clusterContainer, VertexEmbedList *clusterEmbedList, CubicLattice *lattice) :
    ClusterContainer(clusterContainer),
    ClusterEmbedList(clusterEmbedList),
    Lattice(lattice),
    MySubDiagramGenerator(clusterContainer, clusterEmbedList, lattice)
{
    if (this->ClusterContainer->GetL()<2)
        std::cout << "WARNING: Cluster has no subgraphs with two bonds!\n";

    this->FillOneLink(); /// fill list of one-link subdiagrams
    this->AllTotalBondCountsPlusOne[0] = this->OneLink.size()+1; /// set size
    this->LinearIndexMax = this->AllTotalBondCountsPlusOne[0];

    auto tempLambdaHook = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreNNN(index1, index2); }; /// lambda function for hooks
    this->FillBondType(this->Hooks, tempLambdaHook, 3, 2, "HOOK"); /// fill list of hooks
    this->AllTotalBondCountsPlusOne[1] = this->Hooks.size()+1; /// set size
    if (this->Hooks.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[1];

    auto tempLambdaStraightTwoLink = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreThirdNN(index1, index2); }; /// lambda function for straight two-link
    this->FillBondType(this->StraightTwoLink, tempLambdaStraightTwoLink, 3, 2, "STRAIGHT_TWO_LINK");  /// fill list of straight two-link
    this->AllTotalBondCountsPlusOne[2] = this->StraightTwoLink.size()+1; /// set size
    if (this->StraightTwoLink.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[2];

    auto tempLambdaCubeDiagonal = [=](unsigned int index1, unsigned int index2) { return this->Lattice->AreFourthNN(index1, index2); }; /// lambda function for cube diagonal
    this->FillBondType(this->CubeDiagonal, tempLambdaCubeDiagonal, 4, 3, "CUBE_DIAGONAL");  /// fill list of cube diagonal
    this->AllTotalBondCountsPlusOne[3] = this->CubeDiagonal.size()+1; /// set size
    if (this->CubeDiagonal.size()>0) /// add to list of non-zero sizes
        this->LinearIndexMax *= this->AllTotalBondCountsPlusOne[3];

    this->TotalBondCounts = this->OneLink.size()+this->Hooks.size()+this->StraightTwoLink.size()+this->CubeDiagonal.size(); /// precompute total number of bonds (length of valid permutation)
    this->ZCoefficients.resize(this->LinearIndexMax); /// resize
    this->EvaluateZ(); /// evaluate all terms
}

/// get all the subdiagrams corresponding to NN bonds i.e. all \lambda_1 terms
/// ASSUME: cluster has atleast one bond
void ZClusterPureGauge::FillOneLink()
{
    for (int i=0; i<this->MySubDiagramGenerator.GetSizeSubDiagrams(1); ++i)
    {
        VertexEmbedList tempList = this->MySubDiagramGenerator.GetEmbedList(1,i);
#ifdef DEBUG
        if (tempList.GetSize()!=2)
            throw std::invalid_argument("ERROR: FillOneLink accessed a list not of size 2!");
#endif
        std::vector<int> ends(2);
        ends[0] = tempList.begin()->Number;
        ends[1] = std::next(tempList.begin(),1)->Number;
        this->OneLink.push_back(InteractionSubDiagram(ends[0], ends[1], i));
    }
}

/// for a bond type i corresponding to \lambda_i, i>1, identify all subdiagrams corresponding to unique endpoints separated by the appropriate distance
/// @output result: vector onto which we push back identified subdiagrams
/// @param isValidNeighbor: function pointer which takes two cubic lattice indices and returns a true if they are the appropriate type of neighbor (should match up with type of \lambda_i)
/// @param validNbrVertices: number of vertices corresponding to the shortest possible distance between the sites which are neighbors of ith order
/// @param validNbrBonds: taxi-driver distance corresponding to bond of ith order (validNbrVertices=validNbrBonds+1)
/// @param bondType: string denoting bond type (for debugging purposes)
void ZClusterPureGauge::FillBondType(std::vector<InteractionSubDiagram> &result, std::function<bool(unsigned int, unsigned int)> isValidNeighbor, int validNbrVertices, int validNbrBonds, std::string bondType)
{
#ifdef DEBUG
    if (validNbrVertices!=validNbrBonds+1)
        throw std::invalid_argument("ERROR: FillBondType requires validNbrVertices=validNbrBonds+1!\n");
#endif
    if (this->ClusterContainer->GetL()<2)
        return;

    std::set<InteractionSubDiagram> bondSet;
    for (int i=0; i<this->MySubDiagramGenerator.GetSizeSubDiagrams(validNbrBonds); ++i)
    {
        std::vector<int> interactionEndLabels(2,-1);
        if (this->IsSubDiagramValidNhbrPath(this->MySubDiagramGenerator.GetEmbedList(validNbrBonds, i), interactionEndLabels, isValidNeighbor, validNbrVertices))
        {
#ifdef DEBUG
            std::cout << "FillBondType_" << bondType << ": Attempting to add:\n";
            this->MySubDiagramGenerator.GetEmbedList(2, i).PrintList();
#endif
            InteractionSubDiagram tempInteraction(interactionEndLabels[0], interactionEndLabels[1], this->MySubDiagramGenerator.GetSortedLinearIndex(validNbrBonds,i));
            auto result = bondSet.insert(tempInteraction);
#ifdef DEBUG
            if (!result.second)
                std::cout << "FillBondType_" << bondType << ": Equivalent Bond already in list!\n";
            else
                std::cout << "FillBondType_" << bondType << ": New hook found!\n";
#endif
        }
    }
#ifdef DEBUG
    std::cout << "TOTAL_NUMBER_OF_" << bondType << ": " << bondSet.size() << "\n";
    int count = 0;
#endif
    for (auto it=bondSet.begin(); it!=bondSet.end(); ++it)
    {
        result.push_back(*it);
#ifdef DEBUG
        std::cout << bondType << " " << count << " corresponds to subgraph " << it->GetSubDiagramLabel() << " (sorted index) with ends " << it->GetEnd(0) << " and " << it->GetEnd(1) <<  " (labels)\n";
        this->MySubDiagramGenerator.PrintSubDiagram(it->GetSubDiagramLabel());
        count++;
#endif
    }
}

/// for a subgraph of the required taxi-driver distance, see if two of the vertices are ith neighbors! Also return the vertex labels of these "ends"
/// @param list: VertexEmbed list object corresponding to the given subdiagram
/// @output endLabels: vertex labels of the "ends"
/// @param isValidNeighbor: function pointer (are two vertices ith neighbors?)
/// @param validNbrVertices: number of vertices in a taxi-driver path of the appropriate length
bool ZClusterPureGauge::IsSubDiagramValidNhbrPath(const VertexEmbedList& list, std::vector<int>& endLabels, std::function<bool(unsigned int, unsigned int)> isValidNeighbor, int validNbrVertices)
{
#ifdef DEBUG
    if (list.GetSize()!=validNbrVertices)
    {
        std::string errorMessage = "ERROR: IsSubDiagramValidNhbrPath requires list to be of size "+std::to_string(validNbrVertices)+"!\n";
        throw std::invalid_argument(errorMessage);
    }
    if (endLabels.size()!=2)
        throw std::invalid_argument("ERROR: IsSubDiagramValidNhbrPath requires endLabels to be of size 2!\n");
#endif
    for (auto it1=list.begin(); it1!=list.end(); ++it1) /// iterate over pairs of VertexEmbed objects
    {
        for (auto it2=std::next(it1,1); it2!=list.end(); ++it2)
        {
            /// check connectivity
            if (isValidNeighbor(it1->Index, it2->Index)) /// check if they are valid neighbors
            {
#ifdef DEBUG
                std::cout << "IS REQUESTED PATH!\n";
#endif
                /// set labels for ends of hook
                endLabels[0] = it1->Number;
                endLabels[1] = it2->Number;
                return true;
            }
        }
    }
#ifdef DEBUG
    std::cout << "IS NOT REQUESTED PATH!\n";
#endif
    return false;
}

/// evaluate all of the terms in the partition function
/// store results in ZCoefficients
void ZClusterPureGauge::EvaluateZ()
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
        for (int j=0; j<ZClusterPureGauge::NbrCouplings; ++j)
            std::cout << " " << edgesAndBondCount.second[j];
        std::cout << "\n";
        tempContainer.PrintM();
        for (int j=0; j<relabeledEdgesAndVertexMap.second.size(); ++j)
            std::cout << "Vertex " << j+1 << " maps to original vertex " << relabeledEdgesAndVertexMap.second[j] << "\n";
        std::cout << "Weight: " << weight << "\n";
        std::cout << "**** DEBUG_EVALUATEZ ****\n";
#endif
        this->ZCoefficients[this->PowersOfCouplingsToLinearIndex(edgesAndBondCount.second)] += weight;
    }
}

/// print out all non-zero coefficients of Z
void ZClusterPureGauge::PrintZ() const
{
    std::cout << "**** PrintZ ****\n";
    for (int i=0; i<this->LinearIndexMax; ++i)
    {
        double tempCoefficient = this->ZCoefficients[i];
        if (tempCoefficient!=0)
        {
            auto powersCouplings = this->LinearIndexToPowersOfCouplings(i);
            std::cout << "( ";
            for (int j=0; j<ZClusterPureGauge::NbrCouplings; ++j)
                std::cout << powersCouplings[j] << ", ";
            std::cout << "): " << tempCoefficient << "\n";
        }
    }
    std::cout << "**** PrintZ ****\n";
}

/// print the contributions at a given order in the couplings
/// @param powers: requested powers of each coupling \lambda_i
void ZClusterPureGauge::PrintContributionZFixedOrder(const std::array<int, ZClusterPureGauge::NbrCouplings>& powers)
{
    if (!this->ValidPowersOfCouplings(powers))
        throw std::invalid_argument("ERROR: PrintContributionZFixedOrder requires powers to contain integers in the appropriate range!\n");

    std::cout << "**** PrintContributionZFixedOrder ****\n";
    std::cout << "( ";
    for (int j=0; j<ZClusterPureGauge::NbrCouplings; ++j)
        std::cout << powers[j] << ", ";
    std::cout << ")\n";
    for (int i=0; i<this->IntegrandTerms.size(); ++i)
    {
        auto edgesAndBondCount = this->ZIntegrandToUndirectedEdgesAndBondCounts(this->IntegrandTerms[i]); /// get edges and bond counts
        bool correctBondCount = true;
        for (int j=0; j<ZClusterPureGauge::NbrCouplings; ++j)
            if (powers[j]!=edgesAndBondCount.second[j])
                correctBondCount = false;
        if (correctBondCount)
        {
            auto relabeledEdgesAndVertexMap = SubDiagramGenerator::GetRelabeledEdgesAndVertexMap(edgesAndBondCount.first); /// relabeled edges and vertex map
            GraphContainer tempContainer(relabeledEdgesAndVertexMap.second.size(), 1, relabeledEdgesAndVertexMap.first); /// graph container
            PureGaugeWeight tempWeightObject(&tempContainer);
            double weight = tempWeightObject.Weight(); /// calculate weight of graph
            std::cout << "term " << i << " with bond_count:";
            for (int j=0; j<ZClusterPureGauge::NbrCouplings; ++j)
                std::cout << " " << edgesAndBondCount.second[j];
            std::cout << "\n";
            tempContainer.PrintM();
            for (int j=0; j<relabeledEdgesAndVertexMap.second.size(); ++j)
                std::cout << "Vertex " << j+1 << " maps to original vertex " << relabeledEdgesAndVertexMap.second[j] << "\n";
            std::cout << "Weight: " << weight << "\n";
        }
    }
    std::cout << "**** PrintContributionZFixedOrder ****\n";
}

/// convert linear index to a tuple corresponding to the powers of each coupling \lambda_i
/// @param index: linear index
std::array<int, ZClusterPureGauge::NbrCouplings> ZClusterPureGauge::LinearIndexToPowersOfCouplings(int index) const
{
#ifdef DEBUG
    if (index < 0 || index >= this->LinearIndexMax)
        throw std::invalid_argument("ERROR: LinearIndexToPowersOfCouplings requires 0<= index < NTotalBondCounts!\n");
#endif
    std::array<int, ZClusterPureGauge::NbrCouplings> result{0, 0, 0, 0};
    result[0] = index%this->AllTotalBondCountsPlusOne[0]; /// one-link count
    int temp = index/this->AllTotalBondCountsPlusOne[0];
    for (int i=1; i<this->NbrCouplings; ++i)
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
bool ZClusterPureGauge::ValidPowersOfCouplings(const std::array<int, ZClusterPureGauge::NbrCouplings>& powers) const
{
    for (int i=0; i<ZClusterPureGauge::NbrCouplings; ++i)
        if (powers[i] < 0 || powers[i] >= this->AllTotalBondCountsPlusOne[i])
            return false;
    return true;
}

/// convert from a tuple corresponding to powers of couplings to a linear index
/// @param powers: powers of the couplings \lambda_i
int ZClusterPureGauge::PowersOfCouplingsToLinearIndex(const std::array<int, ZClusterPureGauge::NbrCouplings>& powers) const
{
    if (!this->ValidPowersOfCouplings(powers))
        throw std::invalid_argument("ERROR: PowersOfCouplingsToLinearIndex requires powers to contain integers in the appropriate range!\n");

    int linearIndex = powers[0];
    int step = this->AllTotalBondCountsPlusOne[0];
    for (int i=1; i<this->NbrCouplings; ++i)
    {
        if (this->AllTotalBondCountsPlusOne[i]!=0)
        {
            linearIndex += powers[i]*step;
            step *= this->AllTotalBondCountsPlusOne[i];
        }
    }
    return linearIndex;
}

/// recursively generate all combinations of boolean string
/// @param tmp: current string
/// @param nbrBondsRemaining: number of bonds/booleans to add to current string
void ZClusterPureGauge::GenerateTermsIntegrand(std::vector<bool>& tmp, int nbrBondsRemaining)
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

/// convert the boolean string to undirected edges and counts for each bond type
/// @param permutation: boolean strin of length \sum_i N_i
std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGauge::NbrCouplings>> ZClusterPureGauge::ZIntegrandToUndirectedEdgesAndBondCounts(const std::vector<bool>& permutation)
{
#ifdef DEBUG
    if (this->TotalBondCounts!=permutation.size())
        throw std::invalid_argument("ERROR: ZIntegrandToUndirectedEdges requires permutation to be of size TotalBondCounts!\n");
#endif
    std::array<int, ZClusterPureGauge::NbrCouplings> bondCounts{};
    std::vector<UndirectedEdge> result;
    for (int i=0; i<this->OneLink.size(); ++i) /// translate one-link to edges
    {
        if (permutation[i]) /// is this bond requested?
        {
            bondCounts[0]++;
            UndirectedEdge newEdge(this->OneLink[i].GetEnd(0), this->OneLink[i].GetEnd(1));
            result.push_back(newEdge);
        }
    }
    int offset = this->OneLink.size(); /// offset for accessing elements of permutation
    for (int i=0; i<this->Hooks.size(); ++i) /// translate hook to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[1]++;
            UndirectedEdge newEdge(this->Hooks[i].GetEnd(0), this->Hooks[i].GetEnd(1));
            result.push_back(newEdge);
        }
    }
    offset += this->Hooks.size(); /// increment offset
    for (int i=0; i<this->StraightTwoLink.size(); ++i) /// translate straight two-link to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[2]++;
            UndirectedEdge newEdge(this->StraightTwoLink[i].GetEnd(0), this->StraightTwoLink[i].GetEnd(1));
            result.push_back(newEdge);
        }
    }
    offset += this->StraightTwoLink.size(); /// increment offset
    for (int i=0; i<this->CubeDiagonal.size(); ++i) /// translate cube diagonals to edges (endpoints)
    {
        if (permutation[i+offset]) /// is this bond requested?
        {
            bondCounts[3]++;
            UndirectedEdge newEdge(this->CubeDiagonal[i].GetEnd(0), this->CubeDiagonal[i].GetEnd(1));
            result.push_back(newEdge);
        }
    }
#ifdef DEBUG
    for (int i=0; i<ZClusterPureGauge::NbrCouplings; ++i)
        if (bondCounts[i]>this->AllTotalBondCountsPlusOne[i])
            throw std::invalid_argument("ERROR: ZIntegrandToUndirectedEdges requires bondCounts to be less than corresponding element of AllTotalBondCountsPlusOne!\n");
#endif
    return std::pair<std::vector<UndirectedEdge>, std::array<int, ZClusterPureGauge::NbrCouplings>>(result, bondCounts);
}
