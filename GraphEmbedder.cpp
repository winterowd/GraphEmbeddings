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
/// outputFilename: where to output results (two columns: symmetry_factor embedding_number)
void GraphEmbedder::Embed(bool debugJonas)
{

    this->DebugJonas = debugJonas; /// need this member variable set for other member routines

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
    FILE *fpo = fopen(this->Parameters.GetOutputFilename().c_str(), "w");
    if (fpo!=NULL)
        std::cout << "GraphEmbedder::Embed: Successfully opened " << this->Parameters.GetOutputFilename() << std::endl;
    else
        throw std::invalid_argument("GraphEmbedder::Embed: Error opening "+this->Parameters.GetOutputFilename());

    if (this->Parameters.EmbedCorrelator()) /// if we are embedding the two-point function, precompute all lists of initial vertices
    {
        this->GetAllPossiblePairsForFixedVertices(); /// WRONG! FOR NOW EMBED WITH FIRST PAIR IF NOT USING JONAS' INPUT FILE FOR IDENTIFICATION OF ROOTED VERTICES
        //throw std::invalid_argument("Embed: if we embed correlators we need to have debugJonas set!\n");
    }

    std::ifstream myFixedVerticesFile;
    if (this->Parameters.EmbedCorrelator())
    {
        myFixedVerticesFile.open(this->Parameters.GetInputFilenameFixedVertices());
        if (!myFixedVerticesFile.is_open())
            throw std::invalid_argument("GraphEmbedder::Embed: Error opening "+this->Parameters.GetInputFilenameFixedVertices());
    }

    /// read graphs in file
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

        //std::cout << "GraphsFromFile read config " << count << "!\n";
        container.SetGraphFromDenseNauty(g); /// setup container from nauty dense format

        if (this->Parameters.EmbedCorrelator() && this->DebugJonas) /// read fixed vertices from input file
        {
            int v1, v2;
            std::string line;
            if (!std::getline(myFixedVerticesFile, line))
                throw std::invalid_argument("GraphEmbedder::Embed: Error reading "+this->Parameters.GetInputFilenameFixedVertices());
            std::stringstream ss(line);
            if (!(ss >> v1))
                throw std::invalid_argument("GraphEmbedder::Embed: Error parsing first vertex in "+this->Parameters.GetInputFilenameFixedVertices());
            if (!(ss >> v2))
                throw std::invalid_argument("GraphEmbedder::Embed: Error parsing second vertex in "+this->Parameters.GetInputFilenameFixedVertices());
            std::vector<int> tempFixed{v1,v2};
            this->FixedVertexNumbers[0] = tempFixed;
        }

        this->ComputeEmbeddingNumbers(container, fpo);

    }

    DYNFREE(g, g_sz); /// free graph
    DYNFREE(lab, lab_sz); /// free labels
    DYNFREE(ptn, ptn_sz); /// free pattern
    DYNFREE(orbits, orbits_sz); /// free orbits

    fclose(fpo); /// close the file
    myFixedVerticesFile.close(); /// close the file

}

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

/// check if we can nearest neighbor (labeled by nn) to element of list (indexed by elem) is not already occupied
/// list: list of vertices already placed on the lattice
/// elem: index corresponding to already placed vertex that is the nearest neighbor of proposed site
/// nn: index corresponding to NN (lattice dependent)
/// newIndex: (return) lattice index of NN
bool GraphEmbedder::IsProposedNNSiteFree(const std::vector<VertexEmbed> &list, int elem, int nn, int &newIndex)
{
    newIndex = this->Lattice->GetNearestNeighbor(list[elem].Index, nn);
    if (std::find(list.begin(), list.end(), newIndex) != list.end())
        return false;
    return true;
}

/// same as above but with added "degree" of neighbor
bool GraphEmbedder::IsProposedNeighborSiteFree(const VertexEmbedList& list, int elem, int degree, int nn, int &newIndex)
{
    newIndex = this->GetNeighbor(degree, list.GetVertexEmbed(elem).Index, nn);
    if (std::find(list.begin(), list.end(), newIndex) != list.end())
        return false;
    return true;
}

/// TODO: this will need to change when considering edges separated from something other than nearest neighbors
/// check if the proposed site (which we assume to not be occupied) is consistent with vertices already placed on the lattice
/// list: list of vertices already placed on the lattice
/// container: GraphContainer object containing the adjacency matrix
/// elem: index in list corresponding to already placed vertex that is already assumed consistent with proposed site
/// newIndex: lattice index of proposed site
/// newVertexNumber: vertex number that we are placing on the lattice
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

void GraphEmbedder::ComputeEmbeddingNumbers(const GraphContainer& container, FILE *fpo)
{
    this->GetCombinationsOfBonds(container.GetL()); /// generate combos of bonds for the given graph

    for (int i=0; i<this->BondCombinations.size(); ++i)
    {
#ifdef DEBUG
        std::cout << "EMBEDDING WITH COMBO:";
        for (int c=0; c<BondCombinations[i].size(); ++c)
            std::cout << " " << BondCombinations[i][c];
        std::cout << "\n";
#endif
        auto count = this->ComputeEmbeddingNumberCombo(container, this->BondCombinations[i]);

        for (int j=0; j<count.size(); ++j)
            fprintf(fpo, " %d", count[j]); /// print the count of graphs
        for (int d=0; d<this->MaxDegreeNeighbor; ++d)
            fprintf(fpo, " %d", this->BondCombinations[i][d]); /// print bond count
    }
    fprintf(fpo, "\n");
}

/// compute embedding number for a graph (rooted or unrooted)
std::vector<int> GraphEmbedder::ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    auto lists = CreateInitialVertexEmbedLists(container, bondCombo);

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
                                            lists.push_back(temp);  /// add new list to lists

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

        if (lists.size()==0) /// check if anything left?
            break;

    } /// while

    std::vector<int> result;
    if (this->Parameters.EmbedCorrelator())
    {
        std::vector<int> temp(this->MaxDegreeNeighbor, 0);
        for (int i=0; i<lists.size(); ++i)
            temp[lists[i].GetCorrelatorDistanceAsIndex()]++;
        result = temp;
    }
    else
    {

        int temp = 0;
        for (int i=0; i<lists.size(); ++i)
            temp += lists[i].GetNbrChoicesForFirstBond();
        result.push_back(temp);
    }

#ifdef DEBUG
    std::cout << "FOUND " << lists.size() << " embeddings for graph!\n";
    std::cout << "THIS GIVES AN EMBEDDING NUMBER(S) OF:";
    for (int i=0; i<result.size(); ++i)
        std::cout << " " << result[0];
    std::cout << "\n";

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

/// accessors
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

/// got all pairs of vertices to be fixed for two-point function
/// assumes we have N vertices as set in constructor
void GraphEmbedder::GetAllPossiblePairsForFixedVertices()
{
    std::vector<int> tmp;
    std::vector<int> vertices(this->N);
    for (int i = 0; i<vertices.size(); ++i)
        vertices[i] = i+1;
    this->FixedVertexNumbers.clear(); /// clear list of fixed vertices (Note: only need if we call Embed() from main program more than once with different graph orders!)
    this->GenerateUniqueCombinationsWithNoDuplicates(tmp, vertices, 0, 2);
}

/// get all unique combinations of vertices without any duplicate entries using depth first search
void GraphEmbedder::GenerateUniqueCombinationsWithNoDuplicates(std::vector<int>& tmp, const std::vector<int>& vertices, int left, int k)
{
    if (k==0)
    {
        this->FixedVertexNumbers.push_back(tmp);
        return;
    }

    for (int i=left; i<this->N; ++i)
    {
        if (std::find(tmp.begin(), tmp.end(), vertices[i])==tmp.end())
        {
            tmp.push_back(vertices[i]);
            this->GenerateUniqueCombinationsWithNoDuplicates(tmp, vertices, i+1, k-1);
            tmp.pop_back();
        }
    }
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

std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    if (this->Parameters.EmbedCorrelator())
        return CreateInitialVertexEmbedListsRooted(container, bondCombo);
    else
        return CreateInitialVertexEmbedListsNonRooted(container, bondCombo);
}

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

    //this->FixedVertexNumbers.push_back(std::vector<int>{v1,v2}); /// set up fixed vertices

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

    this->FixedVertexNumbers.pop_back(); /// CreateInitialVertexEmbedListsRooted has hardcoded as first element of FixedVertices

}

/// higher-level routine that calls CreateInitialVertexEmbedListsRootedFixed for all pairs of rooted vertices
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRooted(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    /// NOTE: before we have a way to generate weighted graphs we will just take the first set of vertices to be the first ones
    /// DEBUG: Jonas' rooted vertices will be in the FIRST element of FixedVertices
    /// DEBUG: Should only include up to NNN in comparison with Jonas!
    if (this->DebugJonas)
        return this->CreateInitialVertexEmbedListsRootedFixed(container, bondCombo, this->FixedVertexNumbers[0]);

    std::vector<VertexEmbedList> result;
    for (int i=0; i<this->FixedVertexNumbers.size(); ++i) /// loop over vertex pairs
    {
        auto temp = this->CreateInitialVertexEmbedListsRootedFixed(container, bondCombo, this->FixedVertexNumbers[i]);
        for (int j=0; j<temp.size(); ++j)
            result.push_back(temp[j]);
    }
    return result;
}

/// create initial list for each correlator up to MaxDegreeNeighbor
std::vector<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRootedFixed(const GraphContainer& container, const std::vector<int> &bondCombo, const std::vector<int>& fixedVertices)
{
    std::vector<VertexEmbedList> result;

    std::vector<unsigned int> indices(this->Lattice->GetDim(), this->Lattice->GetN()/2); /// spatial indices for first vertex
    auto tempIndex = this->Lattice->GetSiteIndex(indices); /// index for first vertex

    auto v1 = fixedVertices[0]; /// first rooted vertex (color 1)
    auto v2 = fixedVertices[1]; /// second rooted vertex (color 2)
    auto rootedVerticesConnected = container.GetElementAdjacencyMatrix(v1, v2);

    auto d = this->Parameters.GetCorrelatorLength(); /// distance of correlator

    bool bondPossible = true; /// is a bond possible?
    if (d>this->Parameters.GetMaxEmbeddingLength()) /// if correlator distance is greater than max embedding length, no bond is possible
        bondPossible = false;
    else
    {
        if (bondCombo[d]==0) /// if we do not have atleast one step that size in our bond combination, no bond is possible
            bondPossible = false;
    }

    if ((rootedVerticesConnected && bondPossible) || !rootedVerticesConnected) /// we add to list if not vertices NOT connected or if they are connected and we can put a bond between them
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

/// TODO: pass in the bond combination!
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
                (",f", po::value<std::string>(&this->InputFilenameFixedVertices), "input filename containing fixed vertices") /// for debugging with Jonas
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

GraphEmbedderParametersNauty::GraphEmbedderParametersNauty(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}

/// constructor for correlator embedding (rooted graph)
/// takes in what the max embedding length will be and the correlator distance
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxEmbeddingLength, MaxInteractionLength correlatorDistance) :
    NbrChoicesForFirstBond(1),
    TwoPointFunction(true),
    FixedVertices(2, VertexEmbed{-1,-1}),
    CorrelatorDistance(correlatorDistance)
{
#ifdef DEBUG
    std::cout << "Created an rooted graph for a correlator of size " << correlatorDistance << "!\n";
#endif
    switch (maxEmbeddingLength)
    {
    case MaxInteractionLength::NearestNeighbor:
        this->BondCounts.resize(1);
        break;
    case MaxInteractionLength::NextNearestNeighbor:
        this->BondCounts.resize(2);
        break;
    case MaxInteractionLength::ThirdNearestNeighbor:
        this->BondCounts.resize(3);
        break;
    case MaxInteractionLength::FourthNearestNeighbor:
        this->BondCounts.resize(4);
        break;
    default:
        throw std::invalid_argument("Invalid interaction length given to constructor of VertexEmbedList!\n");
        break;
    }
}

/// constructor for unrooted graph
/// takes in what the max interaction length will be
VertexEmbedList::VertexEmbedList(MaxInteractionLength maxLength) :
    TwoPointFunction(false),
    FixedVertices(2, VertexEmbed{-1,-1})
{
#ifdef DEBUG
    std::cout << "Created an unrooted graph!\n";
#endif
    switch (maxLength)
    {
    case MaxInteractionLength::NearestNeighbor:
        this->BondCounts.resize(1);
        break;
    case MaxInteractionLength::NextNearestNeighbor:
        this->BondCounts.resize(2);
        break;
    case MaxInteractionLength::ThirdNearestNeighbor:
        this->BondCounts.resize(3);
        break;
    case MaxInteractionLength::FourthNearestNeighbor:
        this->BondCounts.resize(4);
        break;
    default:
        throw std::invalid_argument("Invalid interaction length given to constructor of VertexEmbedList!\n");
        break;
    }
}

/// add a vertex to the list
void VertexEmbedList::AddVertexEmbed(const VertexEmbed& v)
{
    this->List.push_back(v);
}

/// add a vertex to the list
void VertexEmbedList::AddVertexEmbed(int number, int index)
{
    this->List.push_back(VertexEmbed{number,index});
}

void VertexEmbedList::AddFixedVerticesEmbed(const std::vector<VertexEmbed> &embed)
{
    if (embed.size() != 2)
        throw std::invalid_argument("AddFixedVerticesEmbed requires a vector of size 2!\n");
    if (this->FixedVertices[0].Index!=-1 || this->FixedVertices[1].Index!=-1)
        std::cerr << "WARNING: AddFixedVerticesEmbed fixed vertices already set!\n";
    if (!this->TwoPointFunction)
        std::cerr << "WARNING: AddFixedVerticesEmbed called but TwoPointFunction flag set to false!\n";
    this->FixedVertices = embed;
    this->List.push_back(embed[0]);
    this->List.push_back(embed[1]);
}

/// update the bond count
/// when a vertex is added, we need to check how many times to call this
/// will call atleast once as embedding algorithm will choose a neighbor (NN (0) or NNN (1) etc) of already placed vertex but new vertex could also be adjacent to other already placed vertices
void VertexEmbedList::IncrementBondCount(int dIndex)
{
    if (dIndex < 0 || dIndex >= this->BondCounts.size())
        throw std::invalid_argument("IncrementBondCount requires 0 <= dIndex < BondCounts.size()!\n");
    this->BondCounts[dIndex]++;
}

int VertexEmbedList::GetBondCount(int dIndex) const
{
    if (dIndex < 0 || dIndex >= this->BondCounts.size())
        throw std::invalid_argument("GetBondCount requires 0 <= dIndex < BondCounts.size()!\n");
    return this->BondCounts[dIndex];
}

/// return a vertex in list by reference
VertexEmbed VertexEmbedList::GetVertexEmbed(int index) const
{
    if (index < 0 || index >= this->List.size())
        throw std::invalid_argument("GetVertexEmbed requires 0 <= index < List.size()!\n");
    return this->List[index];
}

/// check if list has repeated vertices
bool VertexEmbedList::HasRepeatedVertices()
{
    for (int i=0; i<(this->List.size()-1); ++i)
        if (std::find_if(std::next(this->List.begin(), i+1), this->List.end(), [this, i](const VertexEmbed& v) { return (this->List[i].Number == v.Number); } ) != this->List.end())
            return true;
    return false;
}

/// check if list has repeated sites
bool VertexEmbedList::HasRepeatedSites()
{
    for (int i=0; i<(this->List.size()-1); ++i)
        if (std::find(std::next(this->List.begin(), i+1), this->List.end(), this->List[i].Index) != this->List.end())
            return true;
    return false;
}

/// accessor for fixed vertices
VertexEmbed VertexEmbedList::GetFixedVertex(int index) const
{
    if (!this->IsTwoPointFunction())
        throw std::logic_error("GetFixedVertex called for unrooted graph!\n");
    if (index !=0 || index != 1)
        throw std::invalid_argument("GetFixedVertex requires index to be 0 or 1!\n");
    return this->FixedVertices[index];
}

/// equality operator for VertexEmbedList
bool operator==(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{

    if (lhs.IsTwoPointFunction()!=rhs.IsTwoPointFunction()) /// check if they are of the same type (rooted or unrooted)
        return false;

    if (lhs.GetSize()!=rhs.GetSize()) /// check that they are the same size
        return false;

    if (lhs.GetNbrBondTypes()!=rhs.GetNbrBondTypes()) /// check that they have the same number of bond types
        return false;

    for (int i=0; i<lhs.GetNbrBondTypes(); ++i) /// check bond counts
        if (lhs.GetBondCount(i)!=rhs.GetBondCount(i))
            return false;

    if (lhs.IsTwoPointFunction()) /// compare fixed vertices if we have a two-point function (order matters!)
        if ((lhs.GetFixedVertex(0) != rhs.GetFixedVertex(0)) || (lhs.GetFixedVertex(1) != rhs.GetFixedVertex(1)))
            return false;

    for (auto it = rhs.begin(); it!=rhs.end(); ++it) /// compare all of the vertices
        if (std::find(lhs.begin(), lhs.end(), *it) == lhs.end())
            return false;

    return true;
}

bool operator!=(const VertexEmbedList& lhs, const VertexEmbedList& rhs)
{
    return !(lhs == rhs);
}

std::ostream& operator<<(std::ostream& os, const VertexEmbedList& list)
{
    os << "********VertexEmbedList********\n";
    for (int i=0; i<list.GetNbrBondTypes(); ++i)
        os << "BondDegree " << i << " with count " << list.GetBondCount(i) << "\n";
    for (int i=0; i<list.GetSize(); ++i)
        os << "Vertex " << list.GetVertexEmbed(i) << "\n";
    return os;
}
