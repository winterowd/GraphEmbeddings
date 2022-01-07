#include "SubDiagramGenerator.h"

/// constructor
/// @param container: pointer to GraphContainer object (managed externally)
/// @param list: pointer to VertexEmbedList object (managed externally)
SubDiagramGenerator::SubDiagramGenerator(GraphContainer *container, VertexEmbedList *list, CubicLattice *lattice) :
    OriginalContainer(container),
    OriginalList(list),
    MyCubicLattice(lattice),
    Vertices(container->GetN()),
    SortedSubDiagramsWithMap(container->GetL(), std::vector<std::pair<int, GraphContainer>>())
{
    /// check that the list has the same size as container
    if (container->GetN()!=list->GetSize())
        throw std::invalid_argument("SubDiagramGenerator requires container and list to be of the same size!\n");

    if (list->IsTwoPointFunction()) /// warn that the subdiagram VertexEmbedLists will not carry info about rooted/unrooted vertices (all unrooted)
        std::cout << "WARNING: SubDiagramGenerator received a two-point funtion! Embeddings of subdiagrams will not contain info about rooted vertices!\n";

    /// construct list of edges
    for (int i=this->OriginalContainer->GetNTimesNMinusOneDiv2()-1; i>=0; --i)
        if (this->OriginalContainer->GetElementAdjacencyMatrix(this->OriginalContainer->GetRowM(i),this->OriginalContainer->GetColM(i)))
            this->Edges.push_back(UndirectedEdge(this->OriginalContainer->GetRowM(i),this->OriginalContainer->GetColM(i)));

    /// list of vertices
    for (int i=0; i<this->OriginalContainer->GetN(); ++i)
        this->Vertices[i] = i+1;

    this->GenerateSubDiagrams(); /// generate subdiagrams
    this->GenerateEmbedListsForSubDiagrams(); /// generate the embed list for the subdiagrams
    this->GenerateCanonicalSubDiagrams();
    this->ComputeDisjointSets(); /// generate the pairwise disjoint sets of connected subgraphs

}

/// generate a list of subdiagrams which are canonical wrt the cubic symmetries and labels
void SubDiagramGenerator::GenerateCanonicalSubDiagrams()
{
    for (int i=0; i<this->NbrSubDiagrams; ++i) /// loop over all subdiagrams
    {
        auto tempCanonical = this->ComputeCanonicalSubgraph(i); /// compute the canonical subgraph (labels and spatial symmetries)
#ifdef DEBUG
        std::cout << "DEBUG_CANONICAL_" << i << "\n";
        std::cout << tempCanonical << "\n";
#endif
        auto tempIt = std::find(this->CanonicalEmbedLists.begin(), this->CanonicalEmbedLists.end(), tempCanonical);
        if (tempIt!=this->CanonicalEmbedLists.end()) /// if it already exists in the list find index and save in vector containing map
        {
            int tempIndex = std::distance(this->CanonicalEmbedLists.begin(), tempIt);
#ifdef DEBUG
            std::cout << "DEBUG_CANONICAL_FOUND element " << tempIndex << "\n";
            std::cout << this->CanonicalEmbedLists[tempIndex] << "\n";
#endif
            this->SubgraphToCanonicalMap.push_back(tempIndex);
        }
        else /// new canonical graph: add to list and update map (last element of current list)
        {
#ifdef DEBUG
            std::cout << "DEBUG_CANONICAL_ADD_NEW!\n";
#endif
            this->CanonicalEmbedLists.push_back(tempCanonical);
            this->SubgraphToCanonicalMap.push_back(this->CanonicalEmbedLists.size()-1); /// corresponds to LAST index  (size-1)
        }
    }
#ifdef DEBUG
    std::cout << "DEBUG_CANONICAL_TOTAL: " << this->CanonicalEmbedLists.size() << "\n";
#endif

}

/// generate the power set recursively
/// Power set of N items made from Power set of N-1 items by taking each set in P_{N-1}, appending the Nth item to it, and appending that new set to P_{N-1}
/// @param elements: set with elements of type T
template<typename T>
std::vector<std::vector<T>> SubDiagramGenerator::GetPowerSet(const std::vector<T>& elements)
{
    if (elements.empty())
        return std::vector<std::vector<T>>(1, std::vector<T>());

    auto allOfThem = GetPowerSet(std::vector<T>(elements.begin()+1,elements.end()));
    T elem = elements[0];

    const int n = allOfThem.size();

    for (int i=0; i<n; ++i)
    {
        const std::vector<T>& s = allOfThem[i];
        allOfThem.push_back(s);
        allOfThem.back().push_back(elem);
    }

    return allOfThem;
}

/// check if the edge is in the given set of vertices i.e. endpoints of edge
/// @param edge: the edge object
/// @param vertices: list of vertices
bool SubDiagramGenerator::IsEdgeInVertexSet(const UndirectedEdge& edge, const std::vector<int>& vertices)
{
#ifdef DEBUG
    if (vertices.size()==0)
    {
        std::cout << "WARNING: In IsEdgeInVertexSet an empty set of vertices given!\n";
        return false;
    }
#endif
    if (vertices.size()<2)
        return false;
    if (std::find(vertices.begin(), vertices.end(), edge.FirstVertex)==vertices.end() || std::find(vertices.begin(), vertices.end(), edge.SecondVertex)==vertices.end())
        return false;
    return true;
}

/// generate all CONNECTED subdiagrams using the power set of edges (assumes original graph is connected)
/// each set of edges inside the power set edges represents a potential subgraph and we test to see if it is connected before adding it
void SubDiagramGenerator::GenerateSubDiagrams()
{
    /// create power set of E
    auto powerE = this->GetPowerSet(this->Edges);

    for (auto eSetIt=(powerE.begin()+1); eSetIt!=powerE.end(); ++eSetIt) /// loop over power sets of E
    {

#ifdef DEBUG
        /// print out set of E
        std::cout << "EdgeSubset: ";
        for (auto elemIt=eSetIt->begin(); elemIt!=eSetIt->end(); ++elemIt)
            std::cout << *elemIt << ", ";
        std::cout << "\n";
#endif
        /// create an edge set that we can send to the constructor of GraphContainer and a valid vertex map
        auto resultRelabelAndMap = this->GetRelabeledEdgesAndVertexMap(*eSetIt);

        GraphContainer subGraph(resultRelabelAndMap.second.size(), 1, resultRelabelAndMap.first);
        if (subGraph.IsConnected())
        {
#ifdef DEBUG
            std::cout << "GRAPH IS CONNECTED!\n";
#endif
            this->AddToSortedSubdiagrams(subGraph, this->VerticesMap.size());
            this->VerticesMap.push_back(resultRelabelAndMap.second);
        }
#ifdef DEBUG
        else
            std::cout << "GRAPH IS DISCONNECTED!\n";
#endif
    }
    int sortedCount = 0;
    for (int i=0; i<this->SortedSubDiagramsWithMap.size(); ++i)
        sortedCount += this->SortedSubDiagramsWithMap[i].size();
    this->NbrSubDiagrams = sortedCount;
#ifdef DEBUG
    std::cout << "NUMBER OF CONNECTED SUBGRAPHS: " << this->NbrSubDiagrams << "\n";
    for (int i=0; i<this->SortedSubDiagramsWithMap.size(); ++i)
        for (int j=0; j<this->SortedSubDiagramsWithMap[i].size(); ++j)
            std::cout << "Diagram of index " << j+1 << " with " << i+1 << " bonds corresponds to unSorted " << this->SortedSubDiagramsWithMap[i][j].first+1 << "\n";
#endif

}

/// compute the canonical subgraph for a given sorted index
/// canonicalize vertex labels with nauty then canonicalize wrt cubic symmetries
/// @param sortedIndex: linear sorted index
CanonicalSubDiagram SubDiagramGenerator::ComputeCanonicalSubgraph(int sortedIndex)
{
#ifdef DEBUG
    std::cout << "Canonicalizing the following subgraph:\n";
    this->PrintSubDiagram(sortedIndex);
#endif
    GraphContainer subgraph = this->GetSubDiagram(sortedIndex);
    auto vertexMap = this->GetVertexMap(sortedIndex);
    VertexEmbedList embedList = this->GetEmbedList(sortedIndex);

    int n = subgraph.GetN();
    int m = SETWORDSNEEDED(n);

    /// compute the canonical graph with NAUTY
    /// declare nauty data structures
    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    /// allocate nauty data structures
    DYNALLOC2(graph, g, g_sz, n, m, "malloc");
    DYNALLOC2(graph, cg, cg_sz, n, m, "malloc");
    DYNALLOC2(int, lab, lab_sz, n, m, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, n, m, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, n, m, "malloc");

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.getcanon = true; /// get canong

    /// set g from container
    subgraph.GetDenseNautyFromGraph(g);

    /// call densenauty
    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, cg);

    /// canonical container
    subgraph.CanonicalRelabeling(lab);

    /// convert vertex labels on embedList (subgraph) from original to canonical relabeled
    VertexEmbedList canonicalList(embedList.GetMaxLength());
    for (auto it=embedList.begin(); it!=embedList.end(); ++it)
    {
        /// find initial relabeling (from original labels to 1,..,N_bsg, where N_bsg is the number of bonds in the subgraph)
        auto tempIt = std::find(vertexMap.begin(), vertexMap.end(), it->Number);
        int tempIndex1 = std::distance(vertexMap.begin(), tempIt);
#ifdef DEBUG
        std::cout << "DEBUG_INVERSE_MAPPING: Vertex " << it->Number << " maps to " << tempIndex1+1 << "\n";
#endif
        /// find canonical relabeling
        int tempIndex2 = -1;
        for (int i=0; i<n; ++i)
        {
            if (lab[i]==tempIndex1)
            {
                tempIndex2 = i;
                break;
            }
        }
        canonicalList.AddVertexEmbed(tempIndex2+1, it->Index); /// add to embed list
#ifdef DEBUG
        std::cout << "DEBUG_CG_MAPPING: Vertex " << tempIndex1+1 << " maps to " << tempIndex2+1 << "\n";
#endif
    }

#ifdef DEBUG
    GraphContainer canonicalSubgraph(n, m, cg);
    if (canonicalSubgraph!=subgraph)
        std::cout << "ERROR: RELABELED_AND_CANONICAL_DIFFER!\n";
#endif

    /// canonicalize with respect to cubic symmetries
    CubicLatticeCanonicalizor canonicalizor(&subgraph, this->MyCubicLattice, canonicalList);

    /// free memory
    DYNFREE(g,g_sz);
    DYNFREE(cg,cg_sz);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    DYNFREE(orbits, orbits_sz);

    return CanonicalSubDiagram(subgraph.GetL(), canonicalizor.GetCanonical());

}

/// create a VertexEmbedList object for each subdiagram using ORIGINAL labels
/// NOTE: EmbedList has same order as VerticesMap
void SubDiagramGenerator::GenerateEmbedListsForSubDiagrams()
{
#ifdef DEBUG
    if (this->NbrSubDiagrams!=this->VerticesMap.size())
        throw std::invalid_argument("ComputeEmbedListForSubDiagram requires VerticesMap to be of size NbrSubDiagrams!\n");
#endif
    for (int i=0; i<this->NbrSubDiagrams; ++i) /// loop over vertex maps
    {
        VertexEmbedList tempList(this->OriginalList->GetMaxLength());
        for (int j=0; j<this->VerticesMap[i].size(); ++j)
            tempList.AddVertexEmbed(this->VerticesMap[i][j], this->GetVertexSiteIndex(this->VerticesMap[i][j]));
        this->EmbedLists.push_back(tempList);
    }
}

/// given a set of edges relabel vertices such that they go from 1,...,V_{sub} where V_{sub} is the number of vertices in the subgraph
/// @param edgeSet: set of edges corresponding to a given subgraph
/// @return std::pair consisting of relabeled set of edges (first) and the vertex map specificing how the new vertex labels map back to the old ones
std::pair<std::vector<UndirectedEdge>, std::vector<int>> SubDiagramGenerator::GetRelabeledEdgesAndVertexMap(const std::vector<UndirectedEdge>& edgeSet)
{
    std::vector<UndirectedEdge> relabeledEdges;
    std::vector<int> vertexMap;
    for (auto it=edgeSet.begin(); it!=edgeSet.end(); ++it) /// loop over edges
    {
        UndirectedEdge tempEdge;

        auto resultFindFirstVertex = std::find(vertexMap.begin(), vertexMap.end(), it->FirstVertex);
        if (resultFindFirstVertex==vertexMap.end()) /// first vertex not yet relabeled (not found in vertexMap)
        {
            vertexMap.push_back(it->FirstVertex); /// v[i] contains the original labeling of vertex label i+1! (NOTE: labels start at 1!)
            tempEdge.FirstVertex = vertexMap.size(); /// relabel first vertex in edge object
#ifdef DEBUG
            std::cout << "Relabeling Vertex " << it->FirstVertex << " as " << vertexMap.size() << "!\n";
#endif
        }
        else /// first vertex already relabeled
        {
#ifdef DEBUG
            std::cout << "Vertex " << it->FirstVertex << " already labeled as " << resultFindFirstVertex-vertexMap.begin()+1 << "!\n";
#endif
            tempEdge.FirstVertex = resultFindFirstVertex-vertexMap.begin()+1; /// relabel first vertex in edge object
        }

        auto resultFindSecondVertex = std::find(vertexMap.begin(), vertexMap.end(), it->SecondVertex);
        if (resultFindSecondVertex==vertexMap.end()) /// second vertex not yet relabeled (not found in vertexMap)
        {
            vertexMap.push_back(it->SecondVertex); /// v[i] contains the original labeling of vertex label i+1! (NOTE: labels start at 1!)
            tempEdge.SecondVertex = vertexMap.size(); /// relabel second vertex in edge object
#ifdef DEBUG
            std::cout << "Relabeling Vertex " << it->SecondVertex << " as " << vertexMap.size() << "!\n";
#endif
        }
        else /// second vertex already relabeled
        {
#ifdef DEBUG
            std::cout << "Vertex " << it->SecondVertex << " already labeled as " << resultFindSecondVertex-vertexMap.begin()+1 << "!\n";
#endif
            tempEdge.SecondVertex = resultFindSecondVertex-vertexMap.begin()+1; /// relabel second vertex in edge
        }
#ifdef DEBUG
        std::cout << "Edge " << *it << " relabeled as " << tempEdge << "\n";
#endif
        relabeledEdges.push_back(tempEdge);
    }
    return std::pair<std::vector<UndirectedEdge>, std::vector<int>>(relabeledEdges, vertexMap);
}

/// get the lattice site index (linear) from a given vertex (using label)
/// @param vertexLabel: vertex for which we want the lattice site index
int SubDiagramGenerator::GetVertexSiteIndex(int vertexLabel) const
{
    if (vertexLabel < 1 || vertexLabel > this->OriginalContainer->GetN())
        throw std::invalid_argument("ERROR: GetVertexSiteIndex requires 1 <= vertexLabel <= N!\n");
    return this->OriginalList->GetVertexSiteIndex(vertexLabel);
}

/// debugging routine
void SubDiagramGenerator::PrintSubDiagrams()
{
    for (int i=0; i<this->NbrSubDiagrams; ++i)
        this->PrintSubDiagram(i);
}

/// print a subdiagram corresponding to the given index (linear SORTED index)
void SubDiagramGenerator::PrintSubDiagram(int index) const
{
    auto convertedIndex = this->IndexConversionSorted(index);
    auto vertexMapIndex = this->SortedSubDiagramsWithMap[convertedIndex.first][convertedIndex.second].first;
    std::cout << "SubDiagram: " << index+1 << "\nVertexMap:\n";
    for (int j=0; j<this->VerticesMap[vertexMapIndex].size(); ++j)
        std::cout << "Vertex " << j+1 << " maps to original vertex " << this->VerticesMap[vertexMapIndex][j] << "\n";
    std::cout << "VertexEmbedList:\n";
    std::cout << this->EmbedLists[vertexMapIndex] << "\n";
    this->SortedSubDiagramsWithMap[convertedIndex.first][convertedIndex.second].second.PrintM();
}

/// see if two subgraphs are disjoint i.e. do not contain a single common vertex (using ORIGINAL labels to search!)
/// @param index1: linear SORTED index for first subdiagram
/// @param index2: linear SORTED index for second subdiagram
bool SubDiagramGenerator::AreDisjoint(int index1, int index2)
{
    if (index1<0 || index1 >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: AreDisjoint requires 0 <= index1 < N_sd!\n");
    if (index2<0 || index2 >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: AreDisjoint requires 0 <= index2 < N_sd!\n");

    /// get indices corresponding to VerticesMap (assume index1 and index2 are in order of SubDiagramsSorted)
    int vertexMapIndex1 = this->GetVertexMapIndexForSubDiagram(index1);
    int vertexMapIndex2 = this->GetVertexMapIndexForSubDiagram(index2);
    auto graphSize = this->GetSubDiagram(index2).GetN();
    for (int i=0; i<graphSize; ++i)
    {
        /// use maps of vertices to see if ith vertex of subdiagram labeled by index2 is also in subdiagram labeled by index1
        if (std::find(this->VerticesMap[vertexMapIndex1].begin(), this->VerticesMap[vertexMapIndex1].end(), this->VerticesMap[vertexMapIndex2][i])!=this->VerticesMap[vertexMapIndex1].end())
            return false;
    }
    return true;
}

/// get size of D_n (number of sets of pairwise disjoint connected subgraphs consisting of n elements)
int SubDiagramGenerator::GetSizeDN(int n) const
{
    if (n<0 || n>=this->NbrSubDiagrams)
        throw std::invalid_argument("GetSizeDN requires 0 <= n < N_sub!\n");
    return this->DisjointSets[n].GetNbrSets();
}

/// compute the D_N iteratively: given D_{N-1}, for each set in D_{N-1} and try to add every subdiagram, checking that all subdiagrams remain pairwise disjoint
void SubDiagramGenerator::ComputeDisjointSets()
{
    SetOfSets<int> D1(1);
    for (int i=0; i<this->NbrSubDiagrams; ++i)
    {
        std::set<int> tempSet;
        tempSet.insert(i);
        D1.AddSet(tempSet);
    }
    this->DisjointSets.push_back(D1);
    int count = 2;
    while (count<=this->NbrSubDiagrams)
    {
        SetOfSets<int> tempDN(count);
        auto dNMinusOne = this->DisjointSets.back(); /// D_{n-1}
#ifdef DEBUG
        std::cout << "D_" << count-1 << " has " << dNMinusOne.GetNbrSets() << " sets!\n";
#endif
        if (dNMinusOne.GetNbrSets()>0) /// empty D_{n-1} means empty D_{n}
        {
            for (int j=0; j<this->NbrSubDiagrams; ++j) /// loop over subdiagrams to add to sets in D_{n-1}
            {
#ifdef DEBUG
                std::cout << "Adding subdiagram " << j+1 << "\n";
                this->PrintSubDiagram(j);
                int setNbr = 0;
#endif
                for (auto itSet=dNMinusOne.begin(); itSet!=dNMinusOne.end(); itSet++)
                {
#ifdef DEBUG
                    std::cout << "adding to set number " << setNbr << " in D_" << count-1 << "\n";
#endif
                    std::set<int> newSet(*itSet); /// save copy of old set
                    bool pairwiseDisjoint = true;
                    for (auto x: *itSet) /// loop over elements in set from D_{n-1}
                    {
                        if (!this->AreDisjoint(x,j)) /// is element disjoint with j?
                        {
                            pairwiseDisjoint = false;
#ifdef DEBUG
                            std::cout << "subdiagram " << j+1 << " not disjoint with element in set number " << setNbr << "!\n";
                            std::cout << "CHECK!\n";
                            this->PrintSubDiagram(x);
#endif
                        }

                    }

                    if (pairwiseDisjoint) /// pairwise disjoint
                    {
#ifdef DEBUG
                        std::cout << "Subdiagram " << j+1 << " is disjoint with all diagrams in set number " << setNbr << " in D_" << count-1 << "\n";
#endif
                        newSet.insert(j); /// add to old set

                        if (newSet.size()!=count)
                            throw std::invalid_argument("ERROR IN ComputeDisjointSets2!\n");

                        tempDN.AddSet(newSet); /// add it to tempDN

                    }
#ifdef DEBUG
                    setNbr++;
#endif
                }
            }
        }
        this->DisjointSets.push_back(tempDN);
        count++; /// increment count
    }
#ifdef DEBUG
    if (this->DisjointSets.size()!=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: DisjointSets should be of size NbrSubDiagrams!\n");
    for (int i=0; i<this->DisjointSets.size(); ++i)
    {
        std::cout << "D_" << i+1 << " has " << this->DisjointSets[i].GetNbrSets() << " sets!\n";
        int count = 0;
        for (auto itSet=this->DisjointSets[i].begin(); itSet!=this->DisjointSets[i].end(); itSet++)
        {
            std::cout << "Subset " << count << " in D_" << i+1 << "!\n";
            for (auto index: *itSet)
                this->PrintSubDiagram(index);
            count++;
        }
    }
#endif
}

/// add a subdiagram along with the unsorted index to SortedSubDiagramsWithMap
/// @param g: GraphContainer object
/// @param indexToUnsorted: index in unsorted vector (VerticesMap, Embedlists)
void SubDiagramGenerator::AddToSortedSubdiagrams(const GraphContainer& g, int indexToUnsorted)
{
    if (g.GetL()>this->OriginalContainer->GetL())
        throw std::invalid_argument("ERROR: AddToSortedSubdiagrams requires g to have number of bonds less than or equal to original graph!\n");
    this->SortedSubDiagramsWithMap[g.GetL()-1].push_back(std::pair<int, GraphContainer>(indexToUnsorted, g));
}

/// accessor for subdiagrams
/// @param index: linear SORTED index
GraphContainer SubDiagramGenerator::GetSubDiagram(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetSubDiagram requires 0 <= index < N_sd!\n");
    auto convertedIndex = this->IndexConversionSorted(sortedIndex);
    return this->SortedSubDiagramsWithMap[convertedIndex.first][convertedIndex.second].second;
}

/// accessor for subdiagrams
/// @param nbrBonds: number of bonds
/// @param graphIndex: index for subdiagrams with given number of bonds
GraphContainer SubDiagramGenerator::GetSubDiagram(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetSubDiagram requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetSubDiagram requires 0 <= graphIndex < N_graphs!\n");
    return this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].second;
}

/// convert from linear SORTED index to (linkIndex, linkSubIndex) where linkIndex = NbrLinks-1 and linkSubIndex is the index within the set of subdiagrams with a given number of links
/// @param index: linear SORTED index
std::pair<int, int> SubDiagramGenerator::IndexConversionSorted(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetVertexMapIndexForSubDiagram requires 0 <= sortedIndex < N_sd!\n");

    int linkIndex = 0;
    for (int i=0; i<this->SortedSubDiagramsWithMap.size(); ++i)
    {
        int temp = sortedIndex-this->SortedSubDiagramsWithMap[i].size(); /// be careful with unsigned integer from size()!
        if (temp>=0)
        {
            linkIndex++;
            sortedIndex = temp;
        }
        else
            break;
    }
    return std::pair<int, int>(linkIndex, sortedIndex);
}

/// get the UNSORTED index given the linear SORTED index!
/// @param index: linear SORTED index
int SubDiagramGenerator::GetVertexMapIndexForSubDiagram(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetVertexMapIndexForSubDiagram requires 0 <= sortedIndex < N_sd!\n");
    auto convertedIndex = this->IndexConversionSorted(sortedIndex);
    return this->SortedSubDiagramsWithMap[convertedIndex.first][convertedIndex.second].first;
}

/// accessor for VerticesMap using the linear SORTED index
/// @param sortedIndex: linear SORTED index
std::vector<int> SubDiagramGenerator::GetVertexMap(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetVertexMap requires 0 <= sortedIndex < N_sd!\n");
    auto unsortedIndex = this->GetVertexMapIndexForSubDiagram(sortedIndex);
    return this->VerticesMap[unsortedIndex];
}

std::vector<int> SubDiagramGenerator::GetVertexMap(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 0 <= graphIndex < N_graphs!\n");
    int unsortedIndex = this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].first;
    return this->VerticesMap[unsortedIndex];
}

/// accessor for EmbedLists using the linear SORTED index
/// @param sortedIndex: linear SORTED index
VertexEmbedList SubDiagramGenerator::GetEmbedList(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetEmbedList requires 0 <= sortedIndex < N_sd!\n");
    auto unsortedIndex = this->GetVertexMapIndexForSubDiagram(sortedIndex);
    return this->EmbedLists[unsortedIndex];
}

VertexEmbedList SubDiagramGenerator::GetEmbedList(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 0 <= graphIndex < N_graphs!\n");
    int unsortedIndex = this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].first;
    return this->EmbedLists[unsortedIndex];
}

int SubDiagramGenerator::GetSortedLinearIndex(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetSortedLinearIndex requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetSortedLinearIndex requires 0 <= graphIndex < N_graphs!\n");
    int result = 0;
    for (int i=1; i<nbrBonds; ++i)
        result += this->GetSizeSubDiagrams(i);
    return result+graphIndex;
}

int SubDiagramGenerator::GetSizeSubDiagrams(int nbrBonds) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetSizeSubDiagrams requires 1 <= nbrBonds <= N_bonds!\n");
    return this->SortedSubDiagramsWithMap[nbrBonds-1].size();
}
