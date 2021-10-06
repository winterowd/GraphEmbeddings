#include "GraphEmbedder.h"

#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// constructor
/// argc: number of command line arguments
/// argv: command line arguments
GraphEmbedder::GraphEmbedder(int argc, char *argv[]) :
    Parameters(argc, argv),
    GetNeighborFunctionPointerArray{&GraphEmbedder::GetNearestNeighbor,&GraphEmbedder::GetNextNearestNeighbor,&GraphEmbedder::GetThirdNearestNeighbor,&GraphEmbedder::GetFourthNearestNeighbor}, /// initialize lists with largest "neighbor" level available
    AreNeighborsFunctionPointerArray{&GraphEmbedder::AreNN,&GraphEmbedder::AreNNN,&GraphEmbedder::AreThirdNN,&GraphEmbedder::AreFourthNN}
{
    switch(this->ResolveLatticeType(this->Parameters.GetLatticeType()))
    {
    case Square:
    {
        this->Lattice = new SquareLattice(100);
        break;
    }
    case Triangular:
    {
        this->Lattice = new TriangularLattice(100);
        break;
    }
    case Cubic:
    {
        this->Lattice = new CubicLattice(100);
        break;
    }
    case Invalid:
    {
        throw std::invalid_argument("GraphEmbedder: requested invalid lattice type "+this->Parameters.GetLatticeType()+"!\n");
    }
    }

    this->N = this->Parameters.GetN();

    for (int v=this->N; v>0; --v)
        this->VertexSet.insert(v);

    this->MWords = SETWORDSNEEDED(this->N); /// set this

    this->MaxDegreeNeighbor = this->Parameters.GetMaxEmbeddingLength()+1; /// set from parameters

    /// store number of neighbors
    this->NbrNeighbors[0] = this->Lattice->GetNbrNN();
    this->NbrNeighbors[1] = this->Lattice->GetNbrNNN();

    nauty_check(WORDSIZE, this->N, this->MWords, NAUTYVERSIONID); /// check if everything is ok

    this->fp = fopen(this->Parameters.GetInputFilename().c_str(), "r");

    if (this->fp!=NULL)
        std::cout << "GraphEmbedder: Successfully opened " << this->Parameters.GetInputFilename() << std::endl;
    else
        throw std::invalid_argument("GraphEmbedder: Error opening "+this->Parameters.GetInputFilename());
}

/// destructor
/// clean up stuff
GraphEmbedder::~GraphEmbedder()
{
    delete this->Lattice; /// delete lattice which we allocated in construtor with new!
    fclose(this->fp); /// close file pointer
}

/// convert string specifying lattice type to private enumeration
GraphEmbedder::LatticeType GraphEmbedder::ResolveLatticeType(const std::string& type)
{
    if (type == "Square") return GraphEmbedder::LatticeType::Square;
    if (type == "Triangular") return GraphEmbedder::LatticeType::Triangular;
    if (type == "Cubic") return GraphEmbedder::LatticeType::Cubic;
    return GraphEmbedder::LatticeType::Invalid;
}

/// wrapper for readg
graph* GraphEmbedder::GetNextGraph(graph *g)
{
    return readg(this->fp,g,0,&this->MWords,&this->N);
}

/// wrapper for densenauty
void GraphEmbedder::CallDenseNauty(graph *g, int *lab, int *ptn, int *orbits, statsblk &stats)
{
    static DEFAULTOPTIONS(options);
    densenauty(g,lab,ptn,orbits,&options,&stats,this->MWords,this->N,NULL);
}

/// public method to call for embedding
/// reads in graphs from a file corresponding to a specific order N and computes the embedding number and size of the automorphism group
void GraphEmbedder::Embed()
{

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc"); /// allocate graph

    /// output file for writing
    FILE *fpo = fopen(this->Parameters.GetOutputFilename().c_str(), "w");
    if (fpo!=NULL)
        std::cout << "GraphEmbedder::Embed: Successfully opened " << this->Parameters.GetOutputFilename() << std::endl;
    else
        throw std::invalid_argument("GraphEmbedder::Embed: Error opening "+this->Parameters.GetOutputFilename());

    /// read graphs in file
    GraphContainer container(this->N, this->MWords);
    int count = 0;
    while (1)
    {
        int symmFactor; /// for rooted graphs
        if (this->Parameters.EmbedCorrelator()) /// correlator?
        {
            char g6temp[MAX_G6_LENGTH];
            int res = fscanf(this->fp, "%s %d\n", g6temp, &symmFactor); /// file format: g6string symm_factor

            if (res == EOF) /// reached the end of file?
                break;
            if (res!=2) /// each line of file should be the same
                throw std::invalid_argument("GraphEmbedder encountered line with wrong format!\n");

            stringtograph(g6temp, g, this->MWords); /// convert string to dense nauty format
        }
        else /// unrooted graph
        {
            graph *gtemp = this->GetNextGraph(g);
            if (gtemp == NULL)
                break;
        }
        count++;
#ifdef DEBUG
        std::cout << "Embed: Read config " << count << "!\n";
#endif
        container.SetGraphFromDenseNauty(g); /// setup container from nauty dense format
        if (this->N != container.GetN()) /// check order!
            std::invalid_argument("Embed found that graph "+std::to_string(count)+" is not of order "+std::to_string(this->N)+"!\n");

        if(this->Parameters.EmbedCorrelator())
            this->ComputeEmbeddingNumbers(container, g, fpo, symmFactor); /// call with symmetry factor passed from file (colored canonicalization)
        else
            this->ComputeEmbeddingNumbers(container, g, fpo); /// compute symmetry factor in routine (uncolored canonicalization)

    }

    DYNFREE(g, g_sz); /// free graph

    fclose(fpo); /// close the file

}

/// debugging routinefor embedding a specific graph (line in file described by graphNbr) for a specific set of bonds (bondCounts)
void GraphEmbedder::EmbedSpecificGraphBondCombo(int graphNbr, const std::vector<int>& bondCounts)
{

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(int,lab,lab_sz); /// declare labels
    DYNALLSTAT(int,ptn,ptn_sz); /// declare pattern
    DYNALLSTAT(int,orbits,orbits_sz); /// declare orbits
    statsblk status; /// status

    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc"); /// allocate graph
    DYNALLOC2(int, lab, lab_sz, this->N, this->MWords, "malloc"); /// allocate labels
    DYNALLOC2(int, ptn, ptn_sz, this->N, this->MWords, "malloc"); /// allocate pattern
    DYNALLOC2(int, orbits, orbits_sz, this->N, this->MWords, "malloc"); /// allocate orbits

    /// output file for writing
    std::string outputFilename = "debug_embed_specific_graph_"+std::to_string(graphNbr)+"_bondCounts";
    for (auto count: bondCounts)
        outputFilename += "_"+std::to_string(count);
    FILE *fpo = fopen(outputFilename.c_str(), "w");
    if (fpo!=NULL)
        std::cout << "GraphEmbedder::EmbedSpecificGraphBondCombo: Successfully opened " << outputFilename << std::endl;
    else
        throw std::invalid_argument("GraphEmbedder::EmbedSpecificGraphBondCombo: Error opening "+outputFilename);

    /// just read first graph for now
    GraphContainer container(this->N, this->MWords);
    int count = 0;
    while (1)
    {
        graph *gtemp = this->GetNextGraph(g);
        if (gtemp == NULL)
            break;
        count++;

        this->CallDenseNauty(g, lab, ptn, orbits, status); /// call densenauty
        writegroupsize(fpo,status.grpsize1,status.grpsize2); /// write group size (symmetry factor) to file

        std::cout << "GraphsFromFile read config " << count << "!\n";
        container.SetGraphFromDenseNauty(g); /// setup container from nauty dense format

        /*if (this->Parameters.GetMaxInteractionLength()==MaxInteractionLength::NearestNeighbor)
        {
            auto embeddingNumber = this->ComputeEmbeddingNumberNN(container); /// TODO: still need to debug this for square lattice (order three graphs)
            fprintf(fpo, " %d\n", embeddingNumber);
        }*/

        if (count == graphNbr)
            this->ComputeEmbeddingNumberCombo(container, bondCounts);

    }

    DYNFREE(g, g_sz); /// free graph
    DYNFREE(lab, lab_sz); /// free labels
    DYNFREE(ptn, ptn_sz); /// free pattern
    DYNFREE(orbits, orbits_sz); /// free orbits

    fclose(fpo); /// close the file
}

/// check if we can add nearest neighbor (labeled by nn) to element of list (indexed by elem) if it is not already occupied
/// @param list: list of vertices already placed on the lattice
/// @param elem: index corresponding to already placed vertex that is the nearest neighbor of proposed site
/// @param nn: index corresponding to NN (lattice dependent)
/// @param newIndex: (return) lattice index of NN
bool GraphEmbedder::IsProposedNNSiteFree(const std::vector<VertexEmbed> &list, int elem, int nn, int &newIndex)
{
    newIndex = this->Lattice->GetNearestNeighbor(list[elem].Index, nn);
    if (std::find(list.begin(), list.end(), newIndex) != list.end())
        return false;
    return true;
}

/// same as above but with added "degree" of neighbor to check for different lengths (NN, NNN, 3N, 4N)
bool GraphEmbedder::IsProposedNeighborSiteFree(const VertexEmbedList& list, int elem, int degree, int nn, int &newIndex)
{
    newIndex = this->GetNeighbor(degree, list.GetVertexEmbed(elem).Index, nn);
    if (std::find(list.begin(), list.end(), newIndex) != list.end())
        return false;
    return true;
}

/// check if the proposed site (which we assume to not be occupied) is consistent with vertices already placed on the lattice
/// @param list: list of vertices already placed on the lattice
/// @param container: GraphContainer object containing the adjacency matrix
/// @param elem: index in list corresponding to already placed vertex that is already assumed consistent with proposed site
/// @param newIndex: lattice index of proposed site
/// @param newVertexNumber: vertex number that we are placing on the lattice
bool GraphEmbedder::IsProposedSiteConsistentWithPreviousVerticesNN(const std::vector<VertexEmbed> &list, const GraphContainer& container, int elem, int newIndex, int newVertexNumber)
{
    for (int i=0; i<list.size(); ++i) /// loop over vertices in set
    {
        if (i != elem) /// assumes we are trying a NN of this vertex (consistent)
        {
            if (container.GetElementAdjacencyMatrix(list[i].Number, newVertexNumber)) /// see if vertex we are trying to place is adjacent to ith vertex in list
            {
                if (!this->Lattice->AreNN(list[i].Index, newIndex)) /// if they are not NN, then it is not consistent
                    return false;
            }
        }
    }
    return true;
}

/// generalization of previous function where bonds can be over any distance allowed by MaxDegreeNeighbor and the given set of bond counts (bondCombo)
bool GraphEmbedder::IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(const VertexEmbedList& list, const GraphContainer &container, const std::vector<int> &bondCombo, int elem, int newIndex, int newVertexNumber, std::vector<int>&  BondCountsToBeAdded)
{
#ifdef DEBUG
    int count = 0;
    for (int d=0; d<this->MaxDegreeNeighbor; ++d) /// check incoming array BondCountsToBeAdded!
    {
        if (BondCountsToBeAdded[d]!=0)
        {
            if (BondCountsToBeAdded[d] < 0 || BondCountsToBeAdded[d] > 1)
                throw std::invalid_argument("IsProposedSiteConsistentWithPreviousVerticesAndBondCounts encountered invalid initial bond count!\n");
            count++;
        }
    }
    if (count > 1)
        throw std::invalid_argument("IsProposedSiteConsistentWithPreviousVerticesAndBondCounts requires BondCountsToBeAdded to have all zeros except one entry which is 1!\n");
#endif

    for (int i=0; i<list.GetSize(); ++i) /// loop over vertices in set
    {
        if (i != elem) /// assumes we are trying a NN of this vertex (consistent)
        {
            auto tempVertexEmbed = list.GetVertexEmbed(i);
            if (container.GetElementAdjacencyMatrix(tempVertexEmbed.Number, newVertexNumber)) /// see if vertex we are trying to place is adjacent to ith vertex in list
            {
                auto allowedBondExists = false; /// two vertices are adjacent but is there an allowed bond on the lattice given the chosen sites?
                for (int d=0; d<this->MaxDegreeNeighbor; ++d) /// loop over "degrees" of neighbors (NN, NNN, etc.)
                {
                    if (this->AreNeighbors(d, tempVertexEmbed.Index, newIndex) && ((list.GetBondCount(d)+BondCountsToBeAdded[d]) < bondCombo[d])) /// are they neighbors of degree d and are there any remaining bonds of length d?
                    {
                        BondCountsToBeAdded[d]++;
                        allowedBondExists = true;
                    }
                }
                if (!allowedBondExists)
                    return false;
            }
        }
    }
    return true;
}


/// Check if lists already contains set of Vertices in toAdd (same vertex numbers and lattice sites)
/// lists: vector of vectors of VertexEmbed objects containing all subgraphs for embedding
/// toAdd: vector of VertexEmbed objects to be added to list
bool GraphEmbedder::IsDuplicate(const std::vector<std::vector<VertexEmbed>>& lists, const std::vector<VertexEmbed>& toAdd)
{
    bool duplicate = false;

    for (int i=0; i<lists.size(); ++i)
    {
        if (lists[i].size() == toAdd.size()) /// only subgraphs with the same number of vertices
        {
            auto ve = std::next(toAdd.begin(),2); /// iterator
            for ( ; ve!=toAdd.end(); ++ve)
            {
#ifdef DEBUG
                std::cout << "Checking element in toAdd against list " << i << "\n";
#endif
                if (std::find(std::next(lists[i].begin(), 2), lists[i].end(), *ve) == lists[i].end()) /// if we cannot find one vertex, lists[i] is not identical to toAdd
                {
#ifdef DEBUG
                    std::cout << "toAdd contains VertexEmbed object with vertex " << ve->Number << " at site " << ve->Index << " which is not in list " << i << "\n";
#endif
                    break;
                }
#ifdef DEBUG
                else
                {
                    std::cout << "toAdd contains VertexEmbed object with vertex " << ve->Number << " at site " << ve->Index << " which is also in list " << i << "\n";
                }
#endif
            }
            if (ve==toAdd.end()) /// if the iterator reached the end, lists[i] is identical to toAdd
            {
                duplicate = true;
#ifdef DEBUG
                std::cout << "toAdd identical to list " << i << "\n";
#endif
                break;
            }
        }
    }
#ifdef DEBUG
    if (!duplicate)
        std::cout << "List can be added! No duplicates!\n";
#endif
    return duplicate;
}

bool GraphEmbedder::AreBondCountsEqual(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{
    if (lhs.GetNbrBondTypes() != rhs.GetNbrBondTypes())
        throw std::invalid_argument("AreBondCountsEqual must compare VertexEmbedList objects with equal nbr of bond types!\n");
    if (lhs.GetNbrBondTypes() != this->MaxDegreeNeighbor)
        throw std::invalid_argument("AreBondCountsEqual only accepts VertexEmbedList objects with nbr of bond types equal to MaxDegreeNeighbor!\n");

    for (int i=0; i<this->MaxDegreeNeighbor; ++i)
        if (lhs.GetBondCount(i)!=rhs.GetBondCount(i))
            return false;
    return true;
}

bool GraphEmbedder::IsDuplicate(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose)
{
    bool result = false;
    for (int i=0; i<lists.size(); ++i)
        if (lists[i]==toAdd)
            result = true;
    return result;
}

bool GraphEmbedder::IsDuplicateOld(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose)
{
    bool duplicate = false;

    for (int i=0; i<lists.size(); ++i)
    {
        /// TODO: also compare bond counts!
        if ((lists[i].GetSize() == toAdd.GetSize()) && this->AreBondCountsEqual(lists[i], toAdd)) /// only subgraphs with the same number of vertices
        {
            auto ve = toAdd.begin(); /// iterator
            for ( ; ve!=toAdd.end(); ++ve)
            {
                if (std::find(lists[i].begin(), lists[i].end(), *ve) == lists[i].end()) /// if we cannot find one vertex, lists[i] is not identical to toAdd
                {

                    break;
                }

            }
            if (ve==toAdd.end()) /// if the iterator reached the end, lists[i] is identical to toAdd
            {
                if (verbose)
                {
                    std::cout << "****** FOUND DUPLICATE ******\n";
                    std::cout << "ALREADY IN LISTS: NUMBER " << i << "\n";
                    std::cout << lists[i] << "\n";
                    std::cout << "TO BE ADDED TO LISTS: \n";
                    std::cout << toAdd << "\n";
                }
                duplicate = true;
                break;
            }
        }
    }
#ifdef DEBUG
    /*if (!duplicate)
        std::cout << "List can be added! No duplicates!\n";
    else
        std::cout << "List cannot be added! Duplicates!\n";*/
    if (duplicate != this->IsDuplicate(lists, toAdd, verbose))
        std::cout << "IsDuplicateNew DOES NOT AGREE with IsDuplicate!\n";
#endif
    return duplicate;
}

/// computes the embedding number of a graph on a given lattice (pointed to by member variable Lattice)
/// container: graph to be embedded
int GraphEmbedder::ComputeEmbeddingNumberNN(const GraphContainer& container)
{
    std::vector<std::vector<VertexEmbed>> lists;

    lists.push_back(this->SelectFirstLinkToEmbedNN(container));

    int vertexCount = 2;

    while (vertexCount < container.GetN())
    {
        auto currentSize = lists.size(); /// add to lists before we remove vectors with the incorrect amount of elements
        for (int i=0; i<currentSize; ++i) /// loop over lists
        {
            auto remainingVertices = this->GetRemainingVertices(lists[i]); /// get vertices remaining to be placed for ith list
#ifdef DEBUG
            if (remainingVertices.size()!=(container.GetN()-vertexCount))
            {
                std::cerr << "Discrepancy with count and number of remaining vertices! " << remainingVertices.size() << " vs " << (container.GetN()-vertexCount) << "\n";
                std::exit(1);
            }
            std::cout << "List " << i << " has yet to add ";
            for (auto it=remainingVertices.begin();  it!=remainingVertices.end(); ++it)
                std::cout << *it << " ";
            std::cout << "\n";
#endif
            for (auto v=remainingVertices.begin(); v!=remainingVertices.end(); ++v) /// loop over vertices to be placed
            {
                for (int j=0; j<lists[i].size(); ++j) /// loop over vertices in ith list
                {
                    /// NOTE: do we just add another loop over degrees? where relevant degree determined for each list in lists (i.e. for each i)
                    for (int nn=1; nn<=this->Lattice->GetNbrNN(); ++nn) /// try placing at nearest neighbor of jth element of ith list
                    {
                        if (container.GetElementAdjacencyMatrix(*v,lists[i][j].Number)) /// is vertex v attached to jth element of ith list?
                        {
                            int nnIndex;
                            if (this->IsProposedNNSiteFree(lists[i], j, nn, nnIndex))
                            {
                                if (this->IsProposedSiteConsistentWithPreviousVerticesNN(lists[i], container, j, nnIndex, *v))
                                {
                                    /// create new list and add to lists
                                    std::vector<VertexEmbed> temp(lists[i]);
                                    temp.push_back(VertexEmbed{*v,nnIndex}); /// add on new vertex to list

                                    /// check if graph already in lists!
                                    if (!this->IsDuplicate(lists, temp))
                                        lists.push_back(temp); /// add new list to lists
                                } /// if consistent
                            } /// if site free
                        } /// if edge exists

                    } /// for nn

                } /// for j
            } /// for v
        } /// for i

        vertexCount++; /// increment vertex count

        /// check that all lists have appropriate number of elements
        lists.erase(std::remove_if(lists.begin(), lists.end(), [&vertexCount](const std::vector<VertexEmbed>& vec) { return vec.size()!=vertexCount; }), lists.end());

        if (lists.size()==0) /// check if anything left?
            break;

    } /// while

#ifdef DEBUG
    std::cout << "FOUND " << lists.size() << " embeddings for graph!\n";

    for (int i=0; i<lists.size(); ++i)
    {
        std::cout << "Embedded graph number " << i << "\n";
        for (int j=0; j<lists[i].size(); ++j)
        {
            std::vector<unsigned int> indices(this->Lattice->GetDim(), 0);
            this->Lattice->GetSiteCoordinates(lists[i][j].Index, indices);
            std::cout << "Vertex " << lists[i][j] << " ";
            for (int k=0; k<this->Lattice->GetDim(); ++k)
                std::cout << indices[k] << " ";
            std::cout << "\n";
        }
    }
#endif

    return this->Lattice->GetNbrNN()*lists.size(); /// multiply by number of NN (first edge fixed by hand)
}

/// test for removing VertexEmbedList objects with wrong size!
void GraphEmbedder::TestEraseWrongSizes(std::vector<VertexEmbedList>& lists, int vertexCount)
{
    lists.erase(std::remove_if(lists.begin(), lists.end(), [&vertexCount](const VertexEmbedList& v) { return v.GetSize()!=vertexCount; }), lists.end());
}

/// compute embedding numbers over all combination of bonds \sum_i N^i_B = NB, where N^i_B is the number of bonds of type i (NN, NNN, 3N 4N) and NB is the total number of bonds
void GraphEmbedder::ComputeEmbeddingNumbers(const GraphContainer& container, graph *g, FILE *fpo, int symmFactor)
{
    if (symmFactor==-1 && this->Parameters.EmbedCorrelator())
        throw std::invalid_argument("ComputeEmbeddingNumbers requires the symmetry factor which is read in from file when embedding correlators (rooted graphs)!\n");

    this->GetCombinationsOfBonds(container.GetL()); /// generate combos of bonds for the given graph

    container.GetDenseNautyFromGraph(g); /// get dense nauty rep of graph
    char *s = ntog6(g,this->MWords,this->N); /// convert to g6

    statsblk status; /// status
    if (!this->Parameters.EmbedCorrelator())
    {
        DYNALLSTAT(int,lab,lab_sz); /// declare labels
        DYNALLSTAT(int,ptn,ptn_sz); /// declare pattern
        DYNALLSTAT(int,orbits,orbits_sz); /// declare orbits

        DYNALLOC2(int, lab, lab_sz, this->N, this->MWords, "malloc"); /// allocate labels
        DYNALLOC2(int, ptn, ptn_sz, this->N, this->MWords, "malloc"); /// allocate pattern
        DYNALLOC2(int, orbits, orbits_sz, this->N, this->MWords, "malloc"); /// allocate orbits

        this->CallDenseNauty(g, lab, ptn, orbits, status); /// call densenauty

        DYNFREE(lab, lab_sz); /// free labels
        DYNFREE(ptn, ptn_sz); /// free pattern
        DYNFREE(orbits, orbits_sz); /// free orbits
    }

    for (int i=0; i<this->BondCombinations.size(); ++i)
    {
#ifdef DEBUG
        std::cout << "EMBEDDING WITH COMBO:";
        for (int c=0; c<BondCombinations[i].size(); ++c)
            std::cout << " " << BondCombinations[i][c];
        std::cout << "\n";
#endif
        auto count = this->ComputeEmbeddingNumberCombo(container, this->BondCombinations[i]);

        fprintf(fpo, "%.*s ", static_cast<int>(strlen(s)-1), s);
        for (int d=0; d<this->MaxDegreeNeighbor; ++d)
            fprintf(fpo, "%d ", this->BondCombinations[i][d]);
        fprintf(fpo, "%d ", count);
        if (!this->Parameters.EmbedCorrelator())
            writegroupsize(fpo,status.grpsize1,status.grpsize2); /// write group size (symmetry factor) to file
        else
            fprintf(fpo, "%d", symmFactor);
        fprintf(fpo, "\n");
    }

}

/// compute embedding number for a graph (rooted or unrooted)
int GraphEmbedder::ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo)
{

    auto lists = CreateInitialVertexEmbedLists(container, bondCombo);

    if (this->Parameters.EmbedCorrelator() && lists.size()==0) /// empty list for correlator i.e. rooted vertices
    {
#ifdef DEBUG
        std::cout << "COULD NOT EMBED ROOTED GRAPH! (ROOTED VERTICES CONNECTED BY A BOND BUT BOND NOT POSSIBLE)\n";
#endif
        return 0; /// return 0 (no point in entering loop as no embedding is possible)
    }

    int vertexCount = 2;

    while (vertexCount < container.GetN())
    {
        auto currentSize = lists.size(); /// add to lists before we remove vectors with the incorrect amount of elements
        //std::cout << "SIZE_AT_TOP_OF_WHILE: " << currentSize << "\n";
        for (int i=0; i<currentSize; ++i) /// loop over lists
        {
            auto remainingVertices = this->GetRemainingVertices(lists[i]); /// get vertices remaining to be placed for ith list
#ifdef DEBUG
            if (remainingVertices.size()!=(container.GetN()-vertexCount))
            {
                std::cerr << "Discrepancy with count and number of remaining vertices! " << remainingVertices.size() << " vs " << (container.GetN()-vertexCount) << "\n";
                std::exit(1);
            }
#endif
            auto dAllowedList = this->GetAllowedBondDegreesOfList(lists[i], bondCombo); /// list of allowed degrees based on bond combo and current VertexEmbedList object
            for (auto v=remainingVertices.begin(); v!=remainingVertices.end(); ++v) /// loop over vertices to be placed
            {
                for (int j=0; j<lists[i].GetSize(); ++j) /// loop over vertices in ith list
                {

                    if (container.GetElementAdjacencyMatrix(*v,lists[i].GetVertexEmbed(j).Number)) /// is vertex v attached to jth element of ith list?
                    {
                        /// loop over allowed degrees of bonds
                        for (int d=0; d<dAllowedList.size(); ++d)
                        {
                            for (int nn=1; nn<=this->GetNbrNeighbors(dAllowedList[d]); ++nn) /// try placing at dth-degree neighbor index nn of jth element of ith list
                            {

                                int nnIndex;
                                if (this->IsProposedNeighborSiteFree(lists[i], j, dAllowedList[d], nn, nnIndex))
                                {

                                    /// adding atleast one bond of degree d (will find out how many others)
                                    std::vector<int> bondCountsToBeAdded(this->MaxDegreeNeighbor,0);
                                    bondCountsToBeAdded[dAllowedList[d]] = 1;
                                    if (this->IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(lists[i], container, bondCombo, j, nnIndex, *v, bondCountsToBeAdded))
                                    {
                                        /// create new list and add to lists
                                        VertexEmbedList temp(lists[i]);
                                        temp.AddVertexEmbed(VertexEmbed{*v,nnIndex}); /// add on new vertex to list

                                        /// update bond counts
                                        this->UpdateBondCounts(temp, bondCountsToBeAdded);

                                        /// check if graph already in lists!
                                        if (!this->IsDuplicate(lists, temp))
                                        {
                                            lists.push_back(temp);  /// add new list to END of lists
#ifdef DEBUG
                                            /// TODO: get coordinates of site corresponding to nnIndex
                                            //std::vector<unsigned int> indices(this->Lattice->GetDim(), 0);
                                            //this->Lattice->GetSiteCoordinates(nnIndex, indices);
                                            //std::cout << "ADDED vertex " << *v << " adjacent to vertex " << lists[i].GetVertexEmbed(j).Number << " at site ";
#endif
                                        }

                                    } /// if consistent
                                } /// if site free

                            } /// for nn
                        } /// for d
                    } /// if edge exists

                } /// for j
            } /// for v

        } /// for i

        vertexCount++; /// increment vertex count

        /// check that all lists have appropriate number of elements
        lists.erase(std::remove_if(lists.begin(), lists.end(), [&vertexCount](const VertexEmbedList& v) { return v.GetSize()!=vertexCount; }), lists.end());

        //std::cout << "ComputeEmbeddingNumberCombo: VertexLists contains " <<  lists.size() << " graphs of size " << vertexCount << "\n";

        if (lists.size()==0) /// check if anything left?
            break;

    } /// while

    int result; /// TODO: why do we return a vector? Shouldn't we return scalar? Do not understand!!!!!!!
    if (this->Parameters.EmbedCorrelator())
    {
        result = lists.size();
    }
    else
    {
        int temp = 0;
        for (int i=0; i<lists.size(); ++i)
            temp += lists[i].GetNbrChoicesForFirstBond();
        result = temp;
    }

#ifdef DEBUG
    std::cout << "FOUND " << lists.size() << " embeddings for graph!\n";
    std::cout << "THIS GIVES AN EMBEDDING NUMBER OF: " << result << "\n";

    for (int i=0; i<lists.size(); ++i)
    {
        std::cout << "Embedded graph number " << i;
        if (this->Parameters.EmbedCorrelator())
            std::cout << " DISTANCE_OF_CORRELATOR " << lists[i].GetCorrelatorDistanceAsIndex();
        std::cout << "\n";
        for (int j=0; j<lists[i].GetSize(); ++j)
        {
            std::vector<unsigned int> indices(this->Lattice->GetDim(), 0);
            this->Lattice->GetSiteCoordinates(lists[i].GetVertexEmbed(j).Index, indices);
            std::cout << "Vertex " << lists[i].GetVertexEmbed(j);
            for (int k=0; k<this->Lattice->GetDim(); ++k)
                std::cout << " " << indices[k];
            std::cout << "\n";
        }
    }
#endif

   return result;
}

void GraphEmbedder::UpdateBondCounts(VertexEmbedList &List, const std::vector<int>& countsToBeAdded)
{
    if (countsToBeAdded.size()!=this->MaxDegreeNeighbor)
        throw std::invalid_argument("UpdateBondCounts requires countsToBeAdded to be of size MaxDegreeNeighbor!\n");
    for (int b=0; b<countsToBeAdded.size(); ++b)
        for (int k=0; k<countsToBeAdded[b]; ++k)
            List.IncrementBondCount(b);
}

/// find first nonzero element in the adjacency matrix (search from the "back") and then place vertices at first NN
/// TODO: this would change if we had anything beyond nearest neighbors: would return a vector of vectors of VertexEmbed objects
std::vector<VertexEmbed> GraphEmbedder::SelectFirstLinkToEmbedNN(const GraphContainer& container)
{
    std::vector<VertexEmbed> result(2); /// return a vector of VertexEmbed objects of length 2

    std::vector<unsigned int> indices(this->Lattice->GetDim(), this->Lattice->GetN()/2); /// spatial indices for first vertex
    auto tempIndex = this->Lattice->GetSiteIndex(indices);
    for (int i=container.GetNTimesNMinusOneDiv2()-1; i>=0; --i)
    {
        if (container.GetElementAdjacencyMatrix(i)) /// first nonzero element of adjacency matrix
        {
            auto v1 = container.GetColM(i); /// find first vertex which corresonds to this edge
            auto v2 = container.GetRowM(i); /// find second vertex which corresonds to this edge
            result[0].Number = v1; result[0].Index = tempIndex;
            result[1].Number = v2; result[1].Index = this->Lattice->GetNearestNeighbor(tempIndex,1); /// second vertex at first NN
            break;
        }
    }
    return result;
}

/// given a vector of VertexEmbed objects, extract the vertex numbers and return an unordered_set containing the vertices that remain
/// listUsedVertices: VertexEmbed objects corresponding to vertices already placed on the graph
std::unordered_set<int> GraphEmbedder::GetRemainingVertices(const std::vector<VertexEmbed>& listUsedVertices)
{
    std::unordered_set<int> usedVertices; /// extract vertices from vector of VertexEmbed objects
    std::for_each(listUsedVertices.begin(), listUsedVertices.end(), [&usedVertices](const VertexEmbed& v) { usedVertices.insert(v.Number); });

    std::unordered_set<int> remainingVertices; /// find vertices that have not been used
    for (std::unordered_set<int>::iterator it=this->VertexSet.begin(); it!=this->VertexSet.end(); ++it)
    {
        auto search = usedVertices.find(*it);
        if (search == usedVertices.end())
            remainingVertices.insert(*it);
    }
    return remainingVertices;
}

/// overloaded version of previous function but with VertexEmbedList object
std::unordered_set<int> GraphEmbedder::GetRemainingVertices(const VertexEmbedList& listUsedVertices)
{
    std::unordered_set<int> usedVertices; /// extract vertices from vector of VertexEmbed objects
    std::for_each(listUsedVertices.begin(), listUsedVertices.end(), [&usedVertices](const VertexEmbed& v) { usedVertices.insert(v.Number); });

    std::unordered_set<int> remainingVertices; /// find vertices that have not been used
    for (std::unordered_set<int>::iterator it=this->VertexSet.begin(); it!=this->VertexSet.end(); ++it)
    {
        auto search = usedVertices.find(*it);
        if (search == usedVertices.end())
            remainingVertices.insert(*it);
    }
    return remainingVertices;
}

std::vector<int> GraphEmbedder::GetAllowedBondDegreesOfList(const VertexEmbedList& list, const std::vector<int> &bondCombo)
{
    std::vector<int> allowed;
    for (int d=0; d<this->MaxDegreeNeighbor; ++d) /// loop through bond degrees
        if (bondCombo[d] > list.GetBondCount(d)) /// if desired bond count is greater than actual bond count for a given degree then we add to allowed counts
            allowed.push_back(d);
    return allowed;
}

/// get the appropriate neighbor
/// degree: level of neighbor (0=NN, 1=NNN, 2=ThirdNN, 3=FourthNN)
/// siteIndex: site index
/// neighborIndex: index of neighbor
int GraphEmbedder::GetNeighbor(int degree, int siteIndex, int neighborIndex)
{
    if (degree < 0 || degree >= NbrLevelsNeighbor)
        throw std::invalid_argument("GetNeighbor requires 0 <= degree < NbrLevelsNeighbor!\n");
    return (this->*GetNeighborFunctionPointerArray[degree])(siteIndex, neighborIndex);
}

/// see if two lattice sites are neighbors
/// degree: level of neighbor (0=NN, 1=NNN, 2=ThirdNN, 3=FourthNN)
/// index1: site index for first site
/// index2: site index for second site
bool GraphEmbedder::AreNeighbors(int degree, int index1, int index2)
{
    if (degree < 0 || degree >= NbrLevelsNeighbor)
        throw std::invalid_argument("AreNeighbors requires 0 <= degree < NbrLevelsNeighbor!\n");
    return (this->*AreNeighborsFunctionPointerArray[degree])(index1, index2);
}

/// accessor
/// degree: level of neighbor (0=NN, 1=NNN; others not yet implemented!)
int GraphEmbedder::GetNbrNeighbors(int degree)
{
    if (degree < 0 || degree >= this->MaxDegreeNeighbor)
        throw std::invalid_argument("GetNbrNeighbors requires 0 <= degree < MaxDegreeNeighbor!\n");
    return this->NbrNeighbors[degree];
}

/// get all combinations of counts for bond types {N_i, i=0(NN),1(NNN),etc. } where \sum_i N_i = nbrBonds
/// nbrBonds: total number of bonds for a particular graph
void GraphEmbedder::GetCombinationsOfBonds(int nbrBonds)
{
    std::vector<int> bondCounts(this->MaxDegreeNeighbor,-1);
    std::vector<int> List(nbrBonds+1);
    for (int i = 0; i<List.size(); ++i)
        List[i] = i;

    this->BondCombinations.clear(); /// clear list of bonds

    this->GenerateCombinations(List, bondCounts, 0); /// generate combinations
}

/// get all combinations of vertices using depth first search
void GraphEmbedder::GenerateCombinations(const std::vector<int>& arr, std::vector<int>& data, int index)
{
    if (index == data.size()) /// created a combo?
    {
        if (std::accumulate(data.begin(), data.end(), 0)==arr.back()) /// combo satisfies \sum_i data_i = arr.back()?
        {
            this->BondCombinations.push_back(data);
        }
        return;
    }

    for (int i=0; i<arr.size(); ++i)
    {
        data[index] = arr[i];
        this->GenerateCombinations(arr, data, index+1);
    }
}

/// calls the appropriate routine for the initial embedded vertices
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    if (this->Parameters.EmbedCorrelator())
        return CreateInitialVertexEmbedListsRooted(container, bondCombo);
    else
        return CreateInitialVertexEmbedListsNonRooted(container, bondCombo);
}

/// debugging routine for initialization of embedding routine for rooted graphs
/// @param filename: file with graph in g6 format
void GraphEmbedder::TestInitialRootedGraphList(std::string filename)
{
    DYNALLSTAT(graph, g, g_sz);
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc");

    fclose(this->fp); /// close main input file
    this->fp = fopen(filename.c_str(), "r"); /// open temp

    if (this->fp!=NULL)
        std::cout << "TestInitialRootedGraphList: Successfully opened " << filename << std::endl;
    else
        throw std::invalid_argument("TestInitialRootedGraphList: Error opening "+filename);

    graph *gtemp = this->GetNextGraph(g);
    if (gtemp == NULL)
    {
        std::cerr << "TestInitialRootedGraphList: Error reading test graph!\n";
        std::exit(1);
    }

    GraphContainer container(this->N, this->MWords, g); /// set up container

    this->GetCombinationsOfBonds(container.GetL()); /// get combos

    int indexCombo = 3; /// hard-code
    if (this->BondCombinations.size()==1)
        indexCombo = 0;

    auto result = this->CreateInitialVertexEmbedListsRooted(container, this->BondCombinations[indexCombo]);

    /// print out result
    /// TODO: each element of result print out

    if (this->fp!=NULL)
        std::cout << "TestInitialRootedGraphList: Successfully opened " << this->Parameters.GetInputFilename() << std::endl;
    else
        throw std::invalid_argument("TestInitialRootedGraphList: Error opening "+this->Parameters.GetInputFilename());

    //this->FixedVertexNumbers.pop_back(); /// CreateInitialVertexEmbedListsRooted has hardcoded as first element of FixedVertices

}

/// creates initial list of embedded vertices for rooted graphs
/// NOTE: embedded vertices ALWAYS taken to be labeled 1 and 2 (0 and 1 in NAUTY convention)! Convention specified in GraphGeneratorNauty
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRooted(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    std::vector<int> rootedVertices{1,2}; /// CONVENTION!
    return this->CreateInitialVertexEmbedListsRootedFixed(container, bondCombo, rootedVertices);
}

/// create initial list for the color given a pair of fixed vertices and a set of bond specifications
/// @param container: graph to be embedded
/// @param bondCombo: set of bond specifications
/// @param rootedVertices: labels of rooted vertices
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRootedFixed(const GraphContainer& container, const std::vector<int> &bondCombo, const std::vector<int>& rootedVertices)
{
    std::vector<VertexEmbedList> result;

    std::vector<unsigned int> indices(this->Lattice->GetDim(), this->Lattice->GetN()/2); /// spatial indices for first vertex
    auto tempIndex = this->Lattice->GetSiteIndex(indices); /// index for first vertex

    auto v1 = rootedVertices[0]; /// first rooted vertex (color 0)
    auto v2 = rootedVertices[1]; /// second rooted vertex (color 1)
    auto rootedVerticesConnected = container.GetElementAdjacencyMatrix(v1, v2);

    auto d = this->Parameters.GetCorrelatorLength(); /// distance of correlator

    bool bondPossible = true; /// is a bond possible?
    if (d>this->Parameters.GetMaxEmbeddingLength()) /// if correlator distance is greater than max embedding length, no bond is possible (NOTE: not sure if this is what we want...)
        bondPossible = false;
    else
    {
        if (bondCombo[d]==0) /// if we do not have atleast one step that size in our bond combination, no bond is possible
            bondPossible = false;
    }

    if (bondPossible || !rootedVerticesConnected) /// we add to list if not vertices NOT connected or if we can place a bond between them
    {
        VertexEmbedList tempList(this->Parameters.GetMaxEmbeddingLength(), static_cast<MaxInteractionLength>(d));
        std::vector<VertexEmbed> tempFixed{VertexEmbed{v1,tempIndex}, VertexEmbed{v2,this->GetNeighbor(d,tempIndex,1)}};
        tempList.AddFixedVerticesEmbed(tempFixed);
        if (rootedVerticesConnected)
            tempList.IncrementBondCount(d);
        result.push_back(tempList);
    }

    return result;
}

/// creates the initial list for non-rooted graphs
/// looks for the first non-zero off-diagonal element of the adjacency matrix and then embeds that bond on the lattice in the canonically defined "first direction" for all allowed bond lengths
/// @param container: graph to be embedded
/// @param bondCombo: set of bond counts
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsNonRooted(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    std::vector<VertexEmbedList> result;

    std::vector<unsigned int> indices(this->Lattice->GetDim(), this->Lattice->GetN()/2); /// spatial indices for first vertex
    auto tempIndex = this->Lattice->GetSiteIndex(indices); /// index for first vertex
    for (int i=container.GetNTimesNMinusOneDiv2()-1; i>=0; --i)
    {
        if (container.GetElementAdjacencyMatrix(i)) /// first nonzero element of adjacency matrix
        {
            /// find vertices which corresonds to this edge
            auto v1 = container.GetColM(i);
            auto v2 = container.GetRowM(i);
            for (int d=0; d<this->MaxDegreeNeighbor; ++d)
            {
                if (bondCombo[d]>0)
                {
                    VertexEmbedList tempList(this->Parameters.GetMaxEmbeddingLength());
                    tempList.AddVertexEmbed(VertexEmbed{v1,tempIndex});
                    tempList.AddVertexEmbed(VertexEmbed{v2,this->GetNeighbor(d,tempIndex,1)});
                    tempList.IncrementBondCount(d);
                    tempList.SetNbrChoicesForFirstBond(this->GetNbrNeighbors(d));
                    result.push_back(tempList);
                }
            }
            break;
        }
    }
    return result;
}

/// process user's command-line arguments
bool GraphEmbedderParametersNauty::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                (",n", po::value<unsigned int>(&this->N)->required(), "Maximum vertex order")
                (",i", po::value<std::string>(&this->InputFilename)->required(), "Input filename")
                (",o", po::value<std::string>(&this->OutputFilename)->required(), "Output filename")
                (",l", po::value<std::string>(&this->LatticeType)->default_value("Cubic"), "Lattice type: Cubic, Square, or Triangular")
                (",d", po::value<MaxInteractionLength>(&this->MaxEmbeddingLength)->default_value(MaxInteractionLength::NearestNeighbor), "MaxEmbeddingLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                (",c", po::value<bool>(&this->Correlator)->default_value(false), "Embed correlators (two-point functions)?")
                (",m", po::value<MaxInteractionLength>(&this->CorrelatorLength)->default_value(MaxInteractionLength::NearestNeighbor), "CorrelatorLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        po::notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    catch (...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }
    return true;
}

/// constructor
GraphEmbedderParametersNauty::GraphEmbedderParametersNauty(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}
