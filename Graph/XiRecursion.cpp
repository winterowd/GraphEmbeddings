#include "XiRecursion.h"

/// constructor
XiRecursionUnrooted::XiRecursionUnrooted(CanonicalGraphManager* manager, const GraphContainer& container, const VertexEmbedList& embedList, CubicLattice *lattice) :
    GraphManager(manager),
    XiContainer(container),
    XiEmbedList(embedList),
    Lattice(lattice)
{
    if (this->GraphManager->GetLMax()<this->XiContainer.GetL())
        throw std::invalid_argument("Error: XiRecursionUnrooted requires that the manager has graphs up to the order of the container for which we compute Xi!\n");

    /// solve the recursion relation for xi
    this->ComputeXiTerms(this->XiContainer, this->XiEmbedList, 1);

    /// clean up elements which have coefficient equal to zero
    this->XiTerms.erase(std::remove_if(this->XiTerms.begin(), this->XiTerms.end(), [](const XiExpansionUnrootedTerm& x) { return (x.GetCoefficient()==0); }), this->XiTerms.end());

    for (int i=0; i<this->XiTerms.size(); ++i)
    {
        std::cout << "FINAL_XI_TERM " << i << "\n";
        std::cout << this->XiTerms[i];
        auto indices = this->XiTerms[i].GetCanonicalCoordinates();
        std::cout << this->GraphManager->GetGraph(indices.first, indices.second); //GetCanonicalSubDiagramContainer(indices.first, indices.second);
    }
}

/// take an embed list which is assumed to be relabeled (vertex labels) subset of original embed list
/// return embed list with original vertex labels (site indices untouched, of course)
/// @param embedList: relabeled embed list
VertexEmbedList XiRecursionUnrooted::RestoreOriginalLabels(const VertexEmbedList& embedList)
{
    if (embedList.GetSize()>this->XiEmbedList.GetSize())
        throw std::invalid_argument("ERROR: RestoreOriginalLabels requires embedList to have a size less than or equal to XiEmbedList!\n");
    VertexEmbedList result(embedList.GetMaxLength());
    for (auto it=embedList.begin(); it!=embedList.end(); ++it)
    {
        auto itOriginal = this->XiEmbedList.begin();
        for (; itOriginal!=this->XiEmbedList.end(); ++itOriginal)
        {
            if (it->Index==itOriginal->Index)
            {
                result.AddVertexEmbed(*itOriginal);
                break;
            }
        }
        if (itOriginal==this->XiEmbedList.end())
            throw std::invalid_argument("ERROR: RestoreOriginalLabels could not find lattice site index in XiEmbedList!\n");
    }
    return result;
}

/// recursive method for computing the terms in the expression
/// \xi(g) = log Z(g) - \sum_{\tilde{g} \in G(g)_{prop sub}} \xi(\tilde{g}), where G_{prop sub} is the set of proper subgraphs of g
/// @param container: graph container
/// @param embedList: embed list for graph
/// @param sign: sign (+/-) which gets distributed in in recursive expression
void XiRecursionUnrooted::ComputeXiTerms(GraphContainer container, VertexEmbedList embedList, int sign)
{
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
        auto temEmbedListFullyCanonical = generator.GetCanonicalSubDiagramEmbedList(1, 0);
        this->AddXiTerm(XiExpansionUnrootedTerm(canonicalIndex, embeddedEdges, temEmbedListFullyCanonical, sign));

        return;
    }

    /// add the graph itself which is a subgraph
    auto tempGraph = generator.GetCanonicalSubDiagramContainer(container.GetL(), 0);
    auto canonicalIndex = this->GraphManager->GetGraphIndex(tempGraph);
    auto edges = generator.GetEmbeddedEdgeSet(container.GetL(), 0);
    auto temEmbedListFullyCanonical = generator.GetCanonicalSubDiagramEmbedList(container.GetL(), 0);
    this->AddXiTerm(XiExpansionUnrootedTerm(canonicalIndex, edges, temEmbedListFullyCanonical, sign));

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
void XiRecursionUnrooted::AddXiTerm(const XiExpansionUnrootedTerm& newTerm)
{
    auto it = std::find(this->XiTerms.begin(), this->XiTerms.end(), newTerm);
    if (it==this->XiTerms.end())
        this->XiTerms.push_back(newTerm);
    else
        it->CombineCoefficients(newTerm.GetCoefficient());
}

/// accessor for Xi terms
XiExpansionUnrootedTerm XiRecursionUnrooted::GetXiTerm(int index) const
{
    if (index<0 || index>=this->XiTerms.size())
        throw std::invalid_argument("ERROR: GetXiTerm requires index to be between 0 and size of XiTerms!\n");
    return this->XiTerms[index];
}
