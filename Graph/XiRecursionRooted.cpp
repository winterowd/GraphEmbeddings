#include "XiRecursionRooted.h"

XiRecursionRooted::XiRecursionRooted(CanonicalGraphManager* manager, const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice *lattice, int embeddingNumber, int maxManhattanDistance) :
    GraphManager(manager),
    XiContainer(container),
    XiEmbedList(embedList),
    Lattice(lattice),
    EmbeddingNumber(embeddingNumber),
    MaxManhattanDistance(maxManhattanDistance)
{
    if (container.GetNbrRooted()!=2)
        throw std::invalid_argument("Error: XiRecursionRooted requires the container to be unrooted!\n");

    if (embedList.GetNbrSetRootedVertices()!=2)
        throw std::invalid_argument("Error: XiRecursionRooted requires the embedList to be unrooted!\n");

    if (this->GraphManager->GetLMax()<this->XiContainer.GetL())
        throw std::invalid_argument("Error: XiRecursionUnrooted requires that the manager has graphs up to the order of the container for which we compute Xi!\n");

    /// solve the recursion relation for xi
    this->ComputeXiTerms(this->XiContainer, this->XiEmbedList, 1);

    /// clean up elements which have coefficient equal to zero
    this->XiTerms.erase(std::remove_if(this->XiTerms.begin(), this->XiTerms.end(), [](const XiExpansionRootedTerm& x) { return (x.GetCoefficient()==0); }), this->XiTerms.end());

    /// compute and store correlators associated with each term in finite-cluster expression
    for (int i=0; i<this->XiTerms.size(); ++i)
    {
        //std::cout << "FINAL_XI_TERM " << i << "\n";
        //std::cout << this->XiTerms[i];
        auto indices = this->XiTerms[i].GetCanonicalCoordinates();
        //std::cout << this->GraphManager->GetRootedGraph(indices.first, indices.second, 2);
        auto tempCorrelator = TwoPointCorrelator(this->GraphManager->GetRootedGraph(indices.first, indices.second, 2), this->XiTerms[i].GetFullyCanonicalEmbedList(), this->Lattice, this->MaxManhattanDistance);
        //tempCorrelator.PrintCorrelatorTerms();
        //std::cout << "DEBUG_EXPAND:\n" << tempCorrelator.GetExpandedCorrelatorGiNaC() << "\n";
        //std::cout << "DEBUG_FULL:\n" << tempCorrelator.GetFullCorrelatorGiNaC() << "\n";
        this->CorrelatorTerms.push_back(tempCorrelator);
    }
}

/// recursive method for computing the terms in the expression
/// \xi(g) = log Z(g) - \sum_{\tilde{g} \in G(g)_{prop sub}} \xi(\tilde{g}), where G_{prop sub} is the set of proper subgraphs of g
/// only keep terms log Z(g) in the expansion if g contains BOTH rooted vertices
/// @param container: graph container
/// @param embedList: embed list for graph
/// @param sign: sign (+/-) which gets distributed in in recursive expression
void XiRecursionRooted::ComputeXiTerms(GraphContainer container, VertexEmbedList embedList, int sign)
{
    if (container.GetNbrRooted()!=embedList.GetNbrSetRootedVertices())
        throw std::logic_error("Error: ComputeXiTerms requires container and embedList to have the same number of rooted vertices!\n");

    /// check if graph is two-rooted; if not, it will not contribute to a correlation function
    if (container.GetNbrRooted()!=2)
        return;

    SubDiagramGenerator generator(container, embedList, this->Lattice);

    if (embedList.GetSize()==2) /// we stop here
    {
        if (generator.GetNbrSubDiagrams()!=1)
            throw std::invalid_argument("ERROR: In ComputeXiTerms base case requires one subdiagram with a single bond!\n");
        if (generator.GetSizeSubDiagrams(1)!=1)
            throw std::invalid_argument("ERROR: In ComputeXiTerms base case requires one subdiagram with a single bond!\n");

        auto tempGraph = generator.GetCanonicalSubDiagramContainer(1, 0);
        auto canonicalIndex = this->GraphManager->GetGraphIndex(tempGraph);
        auto embeddedEdges = generator.GetEmbeddedEdgeSet(1, 0);
        auto tempOriginalList = generator.GetEmbedList(1,0);
        auto tempEmbedListFullyCanonical = generator.GetCanonicalSubDiagramEmbedList(1, 0);
        this->AddXiTerm(XiExpansionRootedTerm(canonicalIndex, embeddedEdges, tempEmbedListFullyCanonical, tempOriginalList, sign));

        return;
    }

    /// add the graph itself which is a subgraph
    auto tempGraph = generator.GetCanonicalSubDiagramContainer(container.GetL(), 0);
    auto canonicalIndex = this->GraphManager->GetGraphIndex(tempGraph);
    auto edges = generator.GetEmbeddedEdgeSet(container.GetL(), 0);
    auto tempEmbedListFullyCanonical = generator.GetCanonicalSubDiagramEmbedList(container.GetL(), 0);
    auto tempOriginalList = generator.GetEmbedList(container.GetL(),0);
    this->AddXiTerm(XiExpansionRootedTerm(canonicalIndex, edges, tempEmbedListFullyCanonical, tempOriginalList, sign));

    for (int l=container.GetL()-1; l>=1; --l) /// loop through all other subgraphs (with bonds less than L) and recurse
    {
        int numGraphsFixedOrder = generator.GetSizeSubDiagrams(l);
        for (int i=0; i<numGraphsFixedOrder; ++i)
            this->ComputeXiTerms(generator.GetCanonicalSubDiagramContainer(l, i), generator.GetEmbedListCanonicalRelabel(l, i), -sign);
    }

    return;
}

/// either add term if we do not already have it or if we do, combine coefficients
/// @param newTerm: new term (log Z(g)) to be added
void XiRecursionRooted::AddXiTerm(const XiExpansionRootedTerm& newTerm)
{
    auto it = std::find(this->XiTerms.begin(), this->XiTerms.end(), newTerm);
    if (it==this->XiTerms.end())
        this->XiTerms.push_back(newTerm);
    else
        it->CombineCoefficients(newTerm.GetCoefficient());
}

/// TODO: write description
XiExpansionRootedTerm XiRecursionRooted::GetXiTerm(int index) const
{
    if (index<0 || index>=this->XiTerms.size())
        throw std::invalid_argument("ERROR: XiRecursionRooted::GetXiTerm index should be between 0 and XiTerms.size()-1!\n");
    return this->XiTerms[index];
}

/// get the full expression for \Xi on the finite-cluster
GiNaC::ex XiRecursionRooted::GetFullXiGiNaC()
{
    GiNaC::ex result;
    for (int i=0; i<this->CorrelatorTerms.size(); ++i)
    {
        auto tempCorrelator = this->CorrelatorTerms[i].GetFullCorrelatorGiNaC();
        //std::cout << "DEBUG_XiRecursionRooted::GetXiGiNaC: " << i << "\n";
        //std::cout << this->XiTerms[i];
        //std::cout << tempCorrelator << "\n";
        result += GiNaC::numeric(this->XiTerms[i].GetCoefficient())*tempCorrelator;
    }
    return result;
}

/// get the expression for \Xi on the finite-cluster expanded to MaxManhattanDistance
GiNaC::ex XiRecursionRooted::GetExpandedXiGiNaC()
{
    GiNaC::ex result;
    for (int i=0; i<this->CorrelatorTerms.size(); ++i)
    {
        auto tempCorrelator = this->CorrelatorTerms[i].GetExpandedCorrelatorGiNaC();
        //std::cout << "DEBUG_XiRecursionRooted::GetXiGiNaC: " << i << "\n";
        //std::cout << this->XiTerms[i];
        //std::cout << tempCorrelator << "\n";
        result += GiNaC::numeric(this->XiTerms[i].GetCoefficient())*tempCorrelator;
    }
    return result;
}

///  compute (embedding/SymmFactor) \times FullGinac
GiNaC::ex XiRecursionRooted::GetFullXiGiNaCWithCoefficient()
{
    if (this->XiContainer.GetSymmFactor()==-1)
        throw std::logic_error("XiRecursionRooted::GetFullXiGiNaCWithCoefficient found XiContainer with unset SymmFactor!\n");
    return GiNaC::numeric(this->EmbeddingNumber,this->XiContainer.GetSymmFactor())*this->GetFullXiGiNaC();
}

/// compute (embedding/SymmFactor) \times ExpandedGinac
GiNaC::ex XiRecursionRooted::GetExpandedXiGiNaCWithCoefficient()
{
    if (this->XiContainer.GetSymmFactor()==-1)
        throw std::logic_error("XiRecursionRooted::GetExpandedXiGiNaCWithCoefficient found XiContainer with unset SymmFactor!\n");
    return GiNaC::numeric(this->EmbeddingNumber,this->XiContainer.GetSymmFactor())*this->GetExpandedXiGiNaC();
}
