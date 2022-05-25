#include "TwoPointCorrelator.h"

/// constructor for two-point correlator object
/// @param container: GraphContainer object for cluster
/// @param embedList: VertexEmbedList corresponding to generalized embedding of graph on CubicLattice (edges of graph DO NOT have to correspond to edges of lattice!)
/// @param lattice: pointer to CubicLattice object
template <typename T>
TwoPointCorrelator<T>::TwoPointCorrelator(const GraphContainer &container, const VertexEmbedList &embedList, CubicLattice* lattice, int maxManhattanDistance, int maxOrderH1, int maxOrderHBar1) :
    ContainerRootedCluster(container),
    EmbedListRootedCluster(embedList),
    Lattice(lattice),
    MaxManhattanDistance(maxManhattanDistance),
    MaxOrderH1(maxOrderH1),
    MaxOrderHBar1(maxOrderHBar1)
{
    /// check number of embedded vertices and check that labels are the same
    if (this->ContainerRootedCluster.GetNbrRooted()!=2)
        throw std::invalid_argument("ERROR: TwoPointCorrelator requires container to be two-rooted!\n");

    if (this->EmbedListRootedCluster.GetNbrSetRootedVertices()!=2)
        throw std::invalid_argument("ERROR: TwoPointCorrelator requires embedList to have both rooted vertices set!\n");

    for (int i=0; i<2; ++i)
    {
        if ((this->ContainerRootedCluster.GetRootedVertex(i)+1)!=this->EmbedListRootedCluster.GetFixedVertex(i).Number)
            throw std::invalid_argument("ERROR: TwoPointCorrelator requires that rooted vertices of container and embedList have same labels!\n");
    }

    /// numerator of connected term at CorrTerms[0]
    this->CorrTerms.push_back(T(this->ContainerRootedCluster, this->EmbedListRootedCluster, this->Lattice, std::vector<bool>{true, false}));

    /// create subgraphs with no rooted graphs and one root (containers and lists)
    GraphContainer containerDenominator(this->ContainerRootedCluster.GetN(), 1, this->ContainerRootedCluster.GetG6String());
    VertexEmbedList embedListDenominator(this->EmbedListRootedCluster.GetMaxLength());
    for (auto it=this->EmbedListRootedCluster.begin(); it!=this->EmbedListRootedCluster.end(); ++it)
        embedListDenominator.AddVertexEmbed(*it);
    /// denominator at CorrTerms[1]
    this->CorrTerms.push_back(T(containerDenominator, embedListDenominator, this->Lattice, std::vector<bool>(0)));

    /// disconnected factor with L at first rooted vertex
    GraphContainer containerDisconnectedL(this->ContainerRootedCluster.GetN(), 1, this->ContainerRootedCluster.GetG6String(), 1);
    containerDisconnectedL.SetRootedVertex(0, this->ContainerRootedCluster.GetRootedVertex(0));
    VertexEmbedList embedListDisconnectedL(this->EmbedListRootedCluster.GetMaxLength(), this->EmbedListRootedCluster.GetCorrelatorDistance());
    embedListDisconnectedL.AddFixedVertexEmbed(0, this->EmbedListRootedCluster.GetFixedVertex(0));
    for (auto it=this->EmbedListRootedCluster.begin();it!=this->EmbedListRootedCluster.end(); ++it)
        embedListDisconnectedL.AddVertexEmbed(*it);
    /// numerator with L at first rooted vertex at CorrTerms[2]
    this->CorrTerms.push_back(T(containerDisconnectedL, embedListDisconnectedL, this->Lattice, std::vector<bool>{true}));

    /// disconnected factor with L* at second rooted vertex
    GraphContainer containerDisconnectedLStar(this->ContainerRootedCluster.GetN(), 1, this->ContainerRootedCluster.GetG6String(), 1);
    containerDisconnectedLStar.SetRootedVertex(0, this->ContainerRootedCluster.GetRootedVertex(0));
    VertexEmbedList embedListDisconnectedLStar(this->EmbedListRootedCluster.GetMaxLength(), this->EmbedListRootedCluster.GetCorrelatorDistance());
    embedListDisconnectedLStar.AddFixedVertexEmbed(0, this->EmbedListRootedCluster.GetFixedVertex(1));
    for (auto it=this->EmbedListRootedCluster.begin();it!=this->EmbedListRootedCluster.end(); ++it)
        embedListDisconnectedLStar.AddVertexEmbed(*it);
    /// numerator with L* at second rooted vertex at CorrTerms[3]
    this->CorrTerms.push_back(T(containerDisconnectedLStar, embedListDisconnectedLStar, this->Lattice, std::vector<bool>{false}));

}

/// create GiNaC::ex object (difference of two rational polynomials in \lambda_i's)
/// CorrTerms[0]/CorrTerms[1] - (CorrTerms[2]/CorrTerms[0])(CorrTerms[3]/CorrTerms[0])
template <typename T>
GiNaC::ex TwoPointCorrelator<T>::GetFullCorrelatorGiNaC()
{
    auto connectedTerm = this->CorrTerms[0].ComputeLambdaPolynomial().GetPolynomial()/this->CorrTerms[1].ComputeLambdaPolynomial().GetPolynomial();
    auto disconnectedTerm = this->CorrTerms[2].ComputeLambdaPolynomial().GetPolynomial()*this->CorrTerms[3].ComputeLambdaPolynomial().GetPolynomial()/(this->CorrTerms[0].ComputeLambdaPolynomial().GetPolynomial()*this->CorrTerms[0].ComputeLambdaPolynomial().GetPolynomial());
    return connectedTerm-disconnectedTerm;
}

/// template specialization for pure gauge
/// create GiNaC::ex object which is the expanded expression for the correlator up to MaxManhattanDistance
template<>
GiNaC::ex TwoPointCorrelator<ZClusterPureGaugeArbEmbedding>::GetExpandedCorrelatorGiNaC()
{
    auto denom = this->CorrTerms[1].ComputeLambdaPolynomial().GetPolynomial();
    auto connectedTerm = AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(this->CorrTerms[0].ComputeLambdaPolynomial().GetPolynomial(), denom, this->MaxManhattanDistance);
    auto disconnectedTerm = AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(this->CorrTerms[2].ComputeLambdaPolynomial().GetPolynomial()*this->CorrTerms[3].ComputeLambdaPolynomial().GetPolynomial(), denom*denom, this->MaxManhattanDistance);
    return connectedTerm-disconnectedTerm;
}

/// template specialization for single flavor static quarks
/// create GiNaC::ex object which is the expanded expression for the correlator up to MaxManhattanDistance
template<>
GiNaC::ex TwoPointCorrelator<ZClusterStaticQuarkArbEmbedding>::GetExpandedCorrelatorGiNaC()
{
    auto denom = this->CorrTerms[1].ComputeLambdaPolynomial().GetPolynomial();
    auto connectedTerm = AuxiliaryRoutinesForGinac::GetFullyExpandedSingleFlavorRationalFunction(this->CorrTerms[0].ComputeLambdaPolynomial().GetPolynomial(), denom, this->MaxManhattanDistance, this->MaxOrderH1, this->MaxOrderHBar1);
    auto disconnectedTerm = AuxiliaryRoutinesForGinac::GetFullyExpandedSingleFlavorRationalFunction(this->CorrTerms[2].ComputeLambdaPolynomial().GetPolynomial()*this->CorrTerms[3].ComputeLambdaPolynomial().GetPolynomial(), denom*denom, this->MaxManhattanDistance, this->MaxOrderH1, this->MaxOrderHBar1);
    return connectedTerm-disconnectedTerm;
}

/// print out each piece of the two-point correlator (debugging purposes)
template <typename T>
void TwoPointCorrelator<T>::PrintCorrelatorTerms()
{
    std::cout << "NUMERATOR_CONNECTED:\n";
    this->CorrTerms[0].PrintZ();
    std::cout << "DENOMINATOR:\n";
    this->CorrTerms[1].PrintZ();
    std::cout << "NUMERATOR_DISCONNECTEDL:\n";
    this->CorrTerms[2].PrintZ();
    std::cout << "NUMERATOR_DISCONNECTEDLStar:\n";
    this->CorrTerms[3].PrintZ();
}

/// enumerate valid template parameters (constructors should have same signature and should have member functions PrintZ and ComputeLambdaPolynomial!)
template class TwoPointCorrelator<ZClusterPureGaugeArbEmbedding>;
template class TwoPointCorrelator<ZClusterStaticQuarkArbEmbedding>;
