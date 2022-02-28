#include "SubDiagramGenerator.h"

/// constructor
/// @param container: pointer to GraphContainer object (managed externally)
/// @param list: pointer to VertexEmbedList object (managed externally)
SubDiagramGenerator::SubDiagramGenerator(const GraphContainer& container, const VertexEmbedList& list, CubicLattice *lattice) :
    OriginalContainer(container),
    OriginalList(list),
    MyCubicLattice(lattice),
    SortedSubDiagramsWithMap(container.GetL(), std::vector<std::pair<int, GraphContainer>>())
{
    /// check that the list has the same size as container
    if (OriginalContainer.GetN()!=OriginalList.GetSize())
        throw std::invalid_argument("SubDiagramGenerator requires container and list to be of the same size!\n");

    if (OriginalList.IsRooted()) /// warn that the subdiagram VertexEmbedLists will not carry info about rooted/unrooted vertices (all unrooted)
        std::cout << "WARNING: SubDiagramGenerator received a two-point funtion! Embeddings of subdiagrams will not contain info about rooted vertices!\n";

     /// check if both container and list are rooted/unrooted
     if (this->OriginalContainer.IsRooted()!=this->OriginalList.IsRooted())
         throw std::invalid_argument("SubDiagramGenerator requires container and list to be both rooted or both unrooted!\n");

    this->GenerateSubDiagrams(); /// generate subdiagrams
    this->GenerateEmbedListsForSubDiagrams(); /// generate the embed list for the subdiagrams
    this->GenerateCanonicalSubDiagrams();

}

/// generate a set of subdiagrams (containers and lists) which are canonical wrt the cubic symmetries and labels
/// each subgraph then maps to an element of the set of canonical subgraphs (i.e. more than one subgraph can map to a single element of the set)
/// NOTE: SubgraphToCanonicalMap is in the same order as
void SubDiagramGenerator::GenerateCanonicalSubDiagrams()
{
    for (int i=0; i<this->NbrSubDiagrams; ++i) /// loop over all subdiagrams
    {
        auto tempCanonicalAndList = this->ComputeCanonicalSubgraphAndList(i); /// compute the canonical subgraph (labels and spatial symmetries)
        this->EmbedListsCanonicalLabels.push_back(tempCanonicalAndList.second); /// add canonically labeled embed list (SORTED ORDER)
#ifdef DEBUG
        std::cout << "DEBUG_CANONICAL_" << i << "\n";
        std::cout << tempCanonicalAndList.first << "\n";
#endif
        auto tempIt = std::find(this->CanonicalSubDiagramList.begin(), this->CanonicalSubDiagramList.end(), tempCanonicalAndList.first);
        if (tempIt!=this->CanonicalSubDiagramList.end()) /// if it already exists in the list find index and save in vector containing map
        {
            int tempIndex = std::distance(this->CanonicalSubDiagramList.begin(), tempIt);
#ifdef DEBUG
            std::cout << "DEBUG_CANONICAL_FOUND element " << tempIndex << "\n";
            std::cout << this->CanonicalSubDiagramList[tempIndex] << "\n";
#endif
            this->SubgraphToCanonicalMap.push_back(tempIndex);
        }
        else /// new canonical graph: add to list and update map (last element of current list)
        {
#ifdef DEBUG
            std::cout << "DEBUG_CANONICAL_ADD_NEW!\n";
            std::cout << tempCanonicalAndList.first << "\n";
#endif
            this->CanonicalSubDiagramList.push_back(tempCanonicalAndList.first);
            this->SubgraphToCanonicalMap.push_back(this->CanonicalSubDiagramList.size()-1); /// corresponds to LAST index  (size-1)
        }
    }
#ifdef DEBUG
    std::cout << "DEBUG_CANONICAL_TOTAL: " << this->CanonicalSubDiagramList.size() << "\n";
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
    auto powerE = this->GetPowerSet(this->OriginalContainer.GetAllEdges());

    for (auto eSetIt=(powerE.begin()+1); eSetIt!=powerE.end(); ++eSetIt) /// loop over power sets of E (skip first element which is the empty set!)
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

        /// create the subgraph (rooted or unrooted)
        auto subGraph = this->CreateSubgraphContainerFromVertexMapAndEdges(resultRelabelAndMap.first, resultRelabelAndMap.second);

        if (subGraph.IsConnected())
        {
#ifdef DEBUG
            std::cout << "GRAPH IS CONNECTED!\n";
            std::cout << "DEBUG_ADDING_GRAPH_TO_SUBDIAGRAMS!\n";
            std::cout << subGraph;
#endif
            this->AddToSortedSubdiagrams(subGraph, this->VerticesMap.size());
            this->VerticesMap.push_back(resultRelabelAndMap.second);
            this->EmbeddedEdgeLists.push_back(this->ConvertUndirectedEdgesToUndirectedEmbeddedEdges(*eSetIt));
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

/// TODO: add comments here
std::pair<CanonicalSubDiagram, VertexEmbedList> SubDiagramGenerator::ComputeCanonicalSubgraphAndList(int sortedIndex)
{
    if (this->GetSubDiagram(sortedIndex).IsRooted())
        return this->ComputeCanonicalSubgraphAndListRooted(sortedIndex);
    else
        return this->ComputeCanonicalSubgraphAndListUnrooted(sortedIndex);
}

/// compute the canonical subgraph for a given sorted index and the map which takes us from canonical vertex labelings to original vertex labelings
/// i.e. vertexMap_cg_to_original[i] contains the original vertex label which is mapped to canonical vertex label (i+1)
/// canonicalize vertex labels with nauty then canonicalize wrt cubic symmetries
/// what we return is a container for the canonical embedded subgraph which
/// @param sortedIndex: linear sorted index for subgraph
std::pair<CanonicalSubDiagram, VertexEmbedList> SubDiagramGenerator::ComputeCanonicalSubgraphAndListUnrooted(int sortedIndex)
{
    GraphContainer subgraph = this->GetSubDiagram(sortedIndex);
    auto vertexMap = this->GetVertexMap(sortedIndex);
    VertexEmbedList embedList = this->GetEmbedList(sortedIndex);

    int n = subgraph.GetN();
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLOC2(int, lab, lab_sz, n, m, "malloc");

    /// get canonical unrooted grpah
    GraphContainer canGraphNew = AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(subgraph.GetN(), subgraph.GetG6String(), lab);

    /// convert vertex labels on embedList (subgraph) from original to canonical relabeled
    VertexEmbedList canonicalList(embedList.GetMaxLength());
    std::vector<int> vertexMapToCG(embedList.GetSize()); /// vertexMapToCG[i] contains original of vertex labeled i+1 after canonicalization
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
        std::cout << "DEBUG_ORIGINAL_TO_CG_MAPPING: Vertex " << it->Number << " maps to " << tempIndex2+1 << "\n";
#endif
        vertexMapToCG[tempIndex2] = it->Number;
    }

    for (int i=0; i<vertexMapToCG.size(); ++i)
        std::cout << "DEBUG_CG_TO_ORIGINAL: VERTEX " << i+1 << " maps to " << vertexMapToCG[i] << "\n";
    std::cout << "DEBUG_CANONICAL_RELABEL_EMBED_LIST:\n";
    std::cout << canonicalList;

    /// canonicalize with respect to cubic symmetries
    CubicLatticeCanonicalizor canonicalizor(canGraphNew, this->MyCubicLattice, canonicalList);

    DYNFREE(lab, lab_sz); /// free label

    return std::pair<CanonicalSubDiagram,VertexEmbedList>(CanonicalSubDiagram(canonicalizor.GetCanonical(), canGraphNew), canonicalList);

}

/// compute the canonical subgraph for a given sorted index and the map which takes us from canonical vertex labelings to original vertex labelings
/// i.e. vertexMap_cg_to_original[i] contains the original vertex label which is mapped to canonical vertex label (i+1)
/// canonicalize vertex labels with nauty then canonicalize wrt cubic symmetries
/// what we return is a container for the canonical embedded subgraph which
/// @param sortedIndex: linear sorted index for subgraph
std::pair<CanonicalSubDiagram, VertexEmbedList> SubDiagramGenerator::ComputeCanonicalSubgraphAndListRooted(int sortedIndex)
{
    GraphContainer subgraph = this->GetSubDiagram(sortedIndex);
    auto vertexMap = this->GetVertexMap(sortedIndex);
    VertexEmbedList embedList = this->GetEmbedList(sortedIndex);

    int n = subgraph.GetN();
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLOC2(int, lab, lab_sz, n, m, "malloc");

    /// canonicalize colored graph (NAUTY)
    GraphContainer canGraphNew = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(subgraph.GetN(), subgraph.GetG6String(), subgraph.GetRootedVertices(), lab);

#ifdef DEBUG
    if (canGraphNew.GetNbrRooted()!=subgraph.GetNbrRooted())
        throw std::logic_error("ERROR: canonical and original graphs have different number of rooted vertices!\n");
#endif

    /// convert vertex labels on embedList (subgraph) from original to canonical relabeled
    VertexEmbedList canonicalList(embedList.GetMaxLength(), embedList.GetCorrelatorDistance());
    std::vector<int> vertexMapToCG(embedList.GetSize()); /// vertexMapToCG[i] contains original of vertex labeled i+1 after canonicalization
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
        if (tempIndex2<canGraphNew.GetNbrRooted()) /// adding a rooted vertex (NAUTY canonicalization so that rooted/colored vertices are relabeled 0 and 1 (1 and 2 if "starting from 1 convention")
            canonicalList.AddFixedVertexEmbed(tempIndex2, tempIndex2+1, it->Index);
        else /// adding an unrooted vertex
            canonicalList.AddVertexEmbed(tempIndex2+1, it->Index); /// add to embed list
#ifdef DEBUG
        std::cout << "DEBUG_CG_MAPPING: Vertex " << tempIndex1+1 << " maps to " << tempIndex2+1 << "\n";
        std::cout << "DEBUG_ORIGINAL_TO_CG_MAPPING: Vertex " << it->Number << " maps to " << tempIndex2+1 << "\n";
#endif
        vertexMapToCG[tempIndex2] = it->Number;
    }

    for (int i=0; i<vertexMapToCG.size(); ++i)
        std::cout << "DEBUG_CG_TO_ORIGINAL: VERTEX " << i+1 << " maps to " << vertexMapToCG[i] << "\n";
    std::cout << "DEBUG_CANONICAL_RELABEL_EMBED_LIST:\n";
    std::cout << canonicalList;

    /// canonicalize with respect to cubic symmetries
    CubicLatticeCanonicalizor canonicalizor(canGraphNew, this->MyCubicLattice, canonicalList);

    DYNFREE(lab, lab_sz); /// free label

    return std::pair<CanonicalSubDiagram,VertexEmbedList>(CanonicalSubDiagram(canonicalizor.GetCanonical(), canGraphNew), canonicalList);

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
        /// get subgraph corresponding to i-th embed list (UNSORTED) to see which vertices are rooted
        auto subGraphContainer = this->GetSubDiagram(this->GetSortedIndexFromUnsortedIndex(i));
        std::cout << "CREATING_EMBEDLIST_FOR_SUBDIAGRAM " << i << "\n";
        if (subGraphContainer.IsRooted())
        {
            VertexEmbedList tempList(this->OriginalList.GetMaxLength(), this->OriginalList.GetCorrelatorDistance());
            /// add rooted vertices
            for (int j=0; j<subGraphContainer.GetNbrRooted(); j++)
            {
                std::cout << "ROOTED_" << j << ": vertex labeled " << subGraphContainer.GetRootedVertex(j)+1 << " which maps to original vertex labeled " << this->VerticesMap[i][subGraphContainer.GetRootedVertex(j)] << "\n";
                tempList.AddFixedVertexEmbed(j, this->VerticesMap[i][subGraphContainer.GetRootedVertex(j)], this->GetVertexSiteIndex(this->VerticesMap[i][subGraphContainer.GetRootedVertex(j)]));
            }
            /// add unrooted
            for (int j=0; j<this->VerticesMap[i].size(); ++j)
            {
                tempList.AddVertexEmbed(this->VerticesMap[i][j], this->GetVertexSiteIndex(this->VerticesMap[i][j]));
            }
            this->EmbedLists.push_back(tempList);
        }
        else
        {
            VertexEmbedList tempList(this->OriginalList.GetMaxLength());
            for (int j=0; j<this->VerticesMap[i].size(); ++j)
                tempList.AddVertexEmbed(this->VerticesMap[i][j], this->GetVertexSiteIndex(this->VerticesMap[i][j]));
            this->EmbedLists.push_back(tempList);
        }
        std::cout << "ADDED_EMBED_LIST: " << i << "\n";
        std::cout << this->EmbedLists[i];
    }
}

/// create container given a set of edges and vertex map
/// takes care of rooted vertices; in subgraph COULD have different color i.e. rooted vertex index 1 could become rooted vertex index 0
/// @param edges: vector of undirected edges
/// @param map: vertex map for relabeling
GraphContainer SubDiagramGenerator::CreateSubgraphContainerFromVertexMapAndEdges(const std::vector<UndirectedEdge>& edges, const std::vector<int>& map)
{
    std::vector<int> rootedNewLabels; /// new labels of rooted vertices (start at 0!)
    for (int i=0; i<this->OriginalContainer.GetNbrRooted(); ++i) /// loop through rooted vertices (0, 1, or 2 rooted vertices)
    {
        auto rootedLabel = this->OriginalContainer.GetRootedVertex(i)+1; /// add one as starts at zero!
        /// find if rootedLabel is in map and if it is, at what index
        auto tempIt = std::find(map.begin(), map.end(), rootedLabel);
        if (tempIt!=map.end())
            rootedNewLabels.push_back(std::distance(map.begin(), tempIt)); /// label from 0,..,N_sb-1
    }
    GraphContainer result(map.size(), 1, edges, rootedNewLabels.size());
    for (int i=0; i<rootedNewLabels.size(); ++i)
        result.SetRootedVertex(i, rootedNewLabels[i]);
    return result;
}

/// given a set of edges relabel vertices such that they go from 1,...,V_{sub} where V_{sub} is the number of vertices in the subgraph
/// @param edgeSet: set of edges corresponding to a given subgraph
/// @return std::pair consisting of relabeled set of edges (first) and the vertex map specificing how the new vertex labels map back to the old ones (second)
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
            tempEdge.FirstVertex = vertexMap.size(); /// relabel first vertex in edge object (starts from 1!)
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
    if (vertexLabel < 1 || vertexLabel > this->OriginalContainer.GetN())
        throw std::invalid_argument("ERROR: GetVertexSiteIndex requires 1 <= vertexLabel <= N!\n");
    return this->OriginalList.GetVertexSiteIndex(vertexLabel);
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
    std::cout << this->EmbedLists[vertexMapIndex];
    std::cout << "Container:\n";
    std::cout << this->SortedSubDiagramsWithMap[convertedIndex.first][convertedIndex.second].second;
    /// TODO: print out canonical container and list
    std::cout << "CanonicalContainer:\n";
    std::cout << this->GetCanonicalSubDiagramContainer(index);
    std::cout << "CanonicalEmbedList:\n";
    std::cout << this->GetCanonicalSubDiagramEmbedList(index);
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
    if (g.GetL()>this->OriginalContainer.GetL())
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

/// accessor for the canonical subgraph corresponding to subgraph at element (nbrBonds, graphIndex)
/// @param nbrBonds: number of bonds
/// @param graphIndex: index for subdiagrams with given number of bonds
GraphContainer SubDiagramGenerator::GetCanonicalSubDiagramContainer(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagram requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagram requires 0 <= graphIndex < N_graphs!\n");
    return this->GetCanonicalSubDiagramContainer(this->GetSortedLinearIndex(nbrBonds, graphIndex)); /// convert to linear index
}

/// accessor for the canonical subgraph corresponding to subgraph at element (sorted linear index)
/// @param sortedIndex: linear sorted index
GraphContainer SubDiagramGenerator::GetCanonicalSubDiagramContainer(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagram requires 0 <= index < N_sd!\n");
    int canonicalIndex = this->SubgraphToCanonicalMap[sortedIndex];
    return this->CanonicalSubDiagramList[canonicalIndex].GetCanonicalContainer();
}

/// accessor for the canonical embed list (wrt NAUTY AND cubic symmetries) corresponding to subgraph at element (nbrBonds, graphIndex)
/// @param nbrBonds: number of bonds
/// @param graphIndex: index for subdiagrams with given number of bonds
VertexEmbedList SubDiagramGenerator::GetCanonicalSubDiagramEmbedList(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagramEmbedList requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagramEmbedList requires 0 <= graphIndex < N_graphs!\n");
    return this->GetCanonicalSubDiagramEmbedList(this->GetSortedLinearIndex(nbrBonds, graphIndex)); /// convert to linear index
}

/// accessor for the canonical embed list (wrt NAUTY AND cubic symmetries) corresponding to subgraph at element (sorted linear index)
/// @param sortedIndex: linear sorted index
VertexEmbedList SubDiagramGenerator::GetCanonicalSubDiagramEmbedList(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetCanonicalSubDiagramEmbedList requires 0 <= index < N_sd!\n");
    int canonicalIndex = this->SubgraphToCanonicalMap[sortedIndex];
    return this->CanonicalSubDiagramList[canonicalIndex].GetCanonicalList();
}

/// convert from linear SORTED index to (linkIndex, linkSubIndex) where linkIndex = NbrLinks-1 and linkSubIndex is the index within the set of subdiagrams with a given number of links
/// @param index: linear SORTED index
std::pair<int, int> SubDiagramGenerator::IndexConversionSorted(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: IndexConversionSorted requires 0 <= sortedIndex < N_sd!\n");

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

/// convert a vector of undirected edges to embedded edges using OriginalList (original VertexEmbedList)
/// @param inputEdges: vector of undirected edges
std::vector<UndirectedEmbeddedEdge> SubDiagramGenerator::ConvertUndirectedEdgesToUndirectedEmbeddedEdges(const std::vector<UndirectedEdge>& inputEdges)
{
    std::vector<UndirectedEmbeddedEdge> result;
    for (auto it=inputEdges.begin(); it!=inputEdges.end(); ++it)
        result.push_back(UndirectedEmbeddedEdge{this->OriginalList.GetVertexSiteIndex(it->FirstVertex), this->OriginalList.GetVertexSiteIndex(it->SecondVertex)});
    return result;
}

/// convert UNSORTED index to linear SORTED index
/// @param unsortedIndex: unsorted index (i.e. embedded lists and vertex maps)
int SubDiagramGenerator::GetSortedIndexFromUnsortedIndex(int unsortedIndex) const
{
    if (unsortedIndex<0 || unsortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetSortedIndexFromUnsortedIndex requires 0 <= sortedIndex < N_sd!\n");
    for (int l=0; l<this->SortedSubDiagramsWithMap.size(); ++l) /// find at what index (l, i) of SortedSubDiagramsWithMap where unsortedIndex is stored
        for (int i=0; i<this->SortedSubDiagramsWithMap[l].size(); ++i)
            if (unsortedIndex==this->SortedSubDiagramsWithMap[l][i].first)
                return this->GetSortedLinearIndex(l+1, i); /// convert (l, i) to linear index (SORTED)
    throw std::invalid_argument("ERROR: In GetSortedIndexFromUnsortedIndex could not convert unsorted index to sorted index!\n");
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

/// accessor for the vertex map between original vertex labels (graph itself) and those of the subgraph
/// @param nbrBonds: number of bonds for the subgraph
/// @param graphIndex: index of subgraph with a particular number of bonds
std::vector<int> SubDiagramGenerator::GetVertexMap(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 0 <= graphIndex < N_graphs!\n");
    int unsortedIndex = this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].first;
    return this->VerticesMap[unsortedIndex];
}

std::vector<UndirectedEmbeddedEdge> SubDiagramGenerator::GetEmbeddedEdgeSet(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetEmbeddedEdgeSet requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetEmbeddedEdgeSet requires 0 <= graphIndex < N_graphs!\n");
    int unsortedIndex = this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].first;
    return this->EmbeddedEdgeLists[unsortedIndex];
}

std::vector<UndirectedEmbeddedEdge> SubDiagramGenerator::GetEmbeddedEdgeSet(int sortedIndex) const
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetEmbeddedEdgeSet requires 0 <= sortedIndex < N_sd!\n");
    auto unsortedIndex = this->GetVertexMapIndexForSubDiagram(sortedIndex);
    return this->EmbeddedEdgeLists[unsortedIndex];
}

VertexEmbedList SubDiagramGenerator::GetEmbedListCanonicalRelabel(int sortedIndex)
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetEmbedListCanonicalRelabel requires 0 <= sortedIndex < N_sd!\n");
    return this->EmbedListsCanonicalLabels[sortedIndex];
}

VertexEmbedList SubDiagramGenerator::GetEmbedListCanonicalRelabel(int nbrBonds, int graphIndex)
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetEmbedListCanonicalRelabel requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetEmbedListCanonicalRelabel requires 0 <= graphIndex < N_graphs!\n");
    return this->GetEmbedListCanonicalRelabel(this->GetSortedLinearIndex(nbrBonds, graphIndex));
}

/// accessor for EmbedLists using the linear SORTED index
/// @param sortedIndex: linear SORTED index
VertexEmbedList SubDiagramGenerator::GetEmbedList(int sortedIndex)
{
    if (sortedIndex<0 || sortedIndex >=this->NbrSubDiagrams)
        throw std::invalid_argument("ERROR: GetEmbedList requires 0 <= sortedIndex < N_sd!\n");
    auto unsortedIndex = this->GetVertexMapIndexForSubDiagram(sortedIndex);
    return this->EmbedLists[unsortedIndex];
}

/// accessor for the embed list of a particular subgraph (vertices retain original labels!)
/// @nbrBonds: number of bonds of a particular subgraph
/// @graphIndex: index of the graph with a particular number of bonds
VertexEmbedList SubDiagramGenerator::GetEmbedList(int nbrBonds, int graphIndex) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 1 <= nbrBonds <= N_bonds!\n");
    if (graphIndex < 0 || graphIndex >= this->SortedSubDiagramsWithMap[nbrBonds-1].size())
        throw std::invalid_argument("ERROR: GetVertexMap requires 0 <= graphIndex < N_graphs!\n");
    int unsortedIndex = this->SortedSubDiagramsWithMap[nbrBonds-1][graphIndex].first;
    return this->EmbedLists[unsortedIndex];
}

/// get the linearized index for the SORTED subdiagram (sorted on number of bonds)
/// @nbrBonds: number of bonds of a particular subgraph
/// @graphIndex: index of the graph with a particular number of bonds
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

/// get the number of subdiagrams with a particular number of bonds
/// @param nbrBonds: number of bonds
int SubDiagramGenerator::GetSizeSubDiagrams(int nbrBonds) const
{
    if (nbrBonds < 1 || nbrBonds > this->SortedSubDiagramsWithMap.size())
        throw std::invalid_argument("ERROR: GetSizeSubDiagrams requires 1 <= nbrBonds <= N_bonds!\n");
    return this->SortedSubDiagramsWithMap[nbrBonds-1].size();
}
