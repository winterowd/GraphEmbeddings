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

    if (!this->Parameters.GetG6().empty())
        this->G6String = this->Parameters.GetG6();

    this->N = this->Parameters.GetN();

    if (this->N!=this->GetGraphSizeFromString(this->G6String))
        throw std::invalid_argument("GraphEmbedder expects N to be equal to the same as g6 graph if given!\n");

    for (int v=this->N; v>0; --v)
        this->VertexSet.insert(v);

    this->MWords = SETWORDSNEEDED(this->N); /// set this

    this->MaxDegreeNeighbor = this->Parameters.GetMaxEmbeddingLength()+1; /// set from parameters

    /// store number of neighbors
    this->NbrNeighbors[0] = this->Lattice->GetNbrNN();
    this->NbrNeighbors[1] = this->Lattice->GetNbrNNN();
    this->NbrNeighbors[2] = this->Lattice->GetNbrThirdNN();
    this->NbrNeighbors[3] = this->Lattice->GetNbrFourthNN();

    nauty_check(WORDSIZE, this->N, this->MWords, NAUTYVERSIONID); /// check if everything is ok

    if (this->Parameters.GetG6().empty()) /// we read in the single g6 string from command line so do not need an input file
    {
        this->fp = fopen(this->Parameters.GetInputFilename().c_str(), "r");

        if (this->fp!=NULL)
            std::cout << "GraphEmbedder: Successfully opened " << this->Parameters.GetInputFilename() << std::endl;
        else
            throw std::invalid_argument("GraphEmbedder: Error opening "+this->Parameters.GetInputFilename());
    }
}

/// destructor
/// clean up stuff
GraphEmbedder::~GraphEmbedder()
{
    delete this->Lattice; /// delete lattice which we allocated in construtor with new!
    if (this->Parameters.GetG6().empty())
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

/// wrapper for densenauty with default options
void GraphEmbedder::CallDenseNauty(graph *g, int *lab, int *ptn, int *orbits, statsblk &stats)
{
    static DEFAULTOPTIONS(options);
    densenauty(g,lab,ptn,orbits,&options,&stats,this->MWords,this->N,NULL);
}

/// public method to call for embedding
void GraphEmbedder::Embed()
{
    if (this->Parameters.GetG6().empty())
        this->EmbedFromFile();
    else
        this->EmbedSingleG6();
}

/// embed a single graph from a g6 string
void GraphEmbedder::EmbedSingleG6()
{

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc"); /// allocate graph

    /// output file for writing
    std::string outputFilename = this->G6String+"_embedding.dat";
    FILE *fpo = fopen(outputFilename.c_str(), "w");
    if (fpo!=NULL)
        std::cout << "GraphEmbedder::Embed: Successfully opened " << outputFilename << std::endl;
    else
        throw std::invalid_argument("GraphEmbedder::Embed: Error opening "+outputFilename);

    this->DenseNautyFromString(this->G6String, g);

    GraphContainer container(this->N, this->MWords, g); /// container from densenauty

    if(this->Parameters.EmbedCorrelator())
    {
        int symmFactor = this->GetSymmFactor(g); /// get symmetry factor!
        this->ComputeEmbeddingNumbers(container, g, fpo, symmFactor);
    }
    else
        this->ComputeEmbeddingNumbers(container, g, fpo);

    DYNFREE(g, g_sz); /// free graph

    fclose(fpo); /// close the file

}

/// reads in graphs from a file corresponding to a specific order N and computes the embedding number and size of the automorphism group
void GraphEmbedder::EmbedFromFile()
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
    //GraphContainer container(this->N, this->MWords);
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
        GraphContainer container(this->N, this->MWords, g); /// setup container from nauty dense format
#ifdef DEBUG
        if (this->N != container.GetN()) /// check order!
            std::invalid_argument("Embed found that graph "+std::to_string(count)+" is not of order "+std::to_string(this->N)+"!\n");
#endif
        if(this->Parameters.EmbedCorrelator())
            this->ComputeEmbeddingNumbers(container, g, fpo, symmFactor); /// call with symmetry factor passed from file (colored canonicalization)
        else
            this->ComputeEmbeddingNumbers(container, g, fpo); /// compute symmetry factor in routine (uncolored canonicalization)

    }

    DYNFREE(g, g_sz); /// free graph

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

/// do the VertexEmbed objects in list contain the index (site) which is the nearest neighbor of some degree of a given vertex
/// @param list: list of vertices which have been embedded
/// @param elem: vertex whose neighboring site we are looking to use to embed a new vertex
/// @param degree: degree of the neighbor (NN, NNN, etc)
/// @param nn: nearest neighbor number
/// @param newIndex (ouput): index of the site at which we will embed the new vertex
bool GraphEmbedder::IsProposedNeighborSiteFree(const VertexEmbedList& list, VertexEmbed elem, int degree, int nn, int &newIndex)
{
    newIndex = this->GetNeighbor(degree, elem.Index, nn);
    if (std::find(list.begin(), list.end(), newIndex) != list.end())
        return false;
    return true;
}

/// check if the proposed NN site (which we assume to not be occupied) is consistent with vertices already placed on the lattice
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
/// @param list: vertices already placed on the lattice
/// @param container: GraphContainer object containing the adjacency matrix
/// @param bondCombo: vector containing the number of each type of bond (NN, NNN, etc)
/// @param elem: vertex in list corresponding to already placed vertex that is already assumed consistent with proposed site
/// @param newIndex: lattice index of proposed site
/// @param newVertexNumber: vertex number that we are placing on the lattice
/// @param BondCountsToBeAdded (output): the new bonds to be counted after placing the vertex on the lattice
bool GraphEmbedder::IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(const VertexEmbedList& list, const GraphContainer &container, const std::vector<int> &bondCombo, VertexEmbed elem, int newIndex, int newVertexNumber, std::vector<int>&  BondCountsToBeAdded)
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

    std::vector<unsigned int> indicesOrigin(this->Lattice->GetDim(), 0);
    auto firstElement = list.begin();
    this->Lattice->GetSiteCoordinates(firstElement->Index, indicesOrigin);

    for (auto it=list.begin(); it!=list.end(); ++it) /// loop over vertices in set
    {
        if (*it != elem) /// assumes we are trying a NN of this vertex (consistent)
        {
            if (container.GetElementAdjacencyMatrix(it->Number, newVertexNumber)) /// see if vertex we are trying to place is adjacent to ith vertex in list
            {
                auto allowedBondExists = false; /// two vertices are adjacent but is there an allowed bond on the lattice given the chosen sites?

                for (int d=0; d<this->MaxDegreeNeighbor; ++d) /// loop over "degrees" of neighbors (NN, NNN, etc.)
                {
                    if (this->AreNeighbors(d, it->Index, newIndex) && ((list.GetBondCount(d)+BondCountsToBeAdded[d]) < bondCombo[d])) /// are they neighbors of degree d and are there any remaining bonds of length d?
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

/// routine which checks if two lists have equal bond counts
/// @param lhs: one list
/// @param rhs: other list
/// @return true if counts are equal, false otherwise
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

/// essentially a wrapper for std::set::find...
bool GraphEmbedder::IsDuplicate(const std::set<VertexEmbedList>& lists, const VertexEmbedList& toAdd)
{
    return (lists.find(toAdd)!=lists.end());
}

/// old std::vector version
bool GraphEmbedder::IsDuplicate(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd)
{
    bool result = false;
    for (int i=0; i<lists.size(); ++i)
        if (lists[i]==toAdd)
            result = true;
    return result;
}

/// oldest version before one could compare VertexEmbedList objects
bool GraphEmbedder::IsDuplicateOld(const std::vector<VertexEmbedList>& lists, const VertexEmbedList& toAdd, bool verbose)
{
    bool duplicate = false;

    for (int i=0; i<lists.size(); ++i)
    {
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
    if (duplicate != this->IsDuplicate(lists, toAdd))
        std::cout << "IsDuplicateNew DOES NOT AGREE with IsDuplicate!\n";
#endif
    return duplicate;
}

/// given a graph G, a list of already embedded vertices and a set of remaining vertices, determine which vertex to embed next and which of embedded vertices it is adjacent to
/// @param container: graph which we are embedding
/// @param embedded: list of embedded vertices
/// @param remainingVertices: set of vertices still to be embedded
/// @return pair consisting of the next vertex to embed and the VertexEmbed object corresponding to its neighbor among the list embedded vertices
std::pair<int, VertexEmbed> GraphEmbedder::DetermineNextVertexToEmbed(const GraphContainer& container, const VertexEmbedList& embedded, const std::unordered_set<int>& remainingVertices)
{
    if (embedded.GetSize()<2)
        throw std::invalid_argument("DetermineNextVertexToEmbed needs argument embedded to be of size greater or equal to 2!\n");
    if (embedded.GetSize()>=container.GetN())
        throw std::invalid_argument("DetermineNextVertexToEmbed needs argument embedded to be of size less than container.GetN()!\n");
    if ((embedded.GetSize()+remainingVertices.size())!=container.GetN())
        throw std::invalid_argument("DetermineNextVertexToEmbed requires container size to equal the sum of the size of embedded and remainingVertices!\n");

    int vNext = -1; /// next vertex to embed
    VertexEmbed vAdjacentNext {-1,-1}; /// already embedded vertex which is adjacent to vNext

    for (auto itRemaining=remainingVertices.begin(); itRemaining!=remainingVertices.end(); ++itRemaining) /// loop over remaining vertices
    {
        for (auto itEmbedded=embedded.begin(); itEmbedded!=embedded.end(); ++itEmbedded)
        {
            if (container.GetElementAdjacencyMatrix(*itRemaining,itEmbedded->Number)) /// found a bond between one of the remaining vertices and one of the already embedded vertices!
            {
                vNext = *itRemaining;
                vAdjacentNext = *itEmbedded;
                goto endOfLoop; /// break out of nested for loop
            }
        }
    }
    endOfLoop:
    if (vNext == -1)
        throw std::invalid_argument("GraphEmbedder::DetermineNextVertexToEmbed did not find a remaining vertex which is connected to previous vertices!\n");
#ifdef DEBUG
    std::cout << "Inside DetermineNextVertexToEmbed\n";
    this->PrintVertexEmbedList(embedded);
    std::cout << "Going to embed vertex " << vNext << " what is adjacent to already embedded vertex " << vAdjacentNext.Number << "!\n";
#endif
    return std::pair<int,VertexEmbed>(vNext, vAdjacentNext);
}

/// old version which computes the embedding number of a graph on a given lattice (NN only)
/// @param container: graph to be embedded
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

    this->GetCombinationsOfBondsFixedNumberOfBonds(container.GetL()); /// generate combos of bonds for fixed bond number
    /// TODO: how to incorporate the constraint with Manhattan distance? think about this...

    container.GetDenseNautyFromGraph(g); /// get dense nauty rep of graph
    char *s = ntog6(g,this->MWords,this->N); /// convert to g6

    if (symmFactor!=-1) /// get the symmetry factor if default value
        symmFactor = this->GetSymmFactor(g);

    for (int i=0; i<this->BondCombinations.size(); ++i)
    {
#ifdef DEBUG
        std::cout << "EMBEDDING WITH COMBO:";
        for (int c=0; c<BondCombinations[i].size(); ++c)
            std::cout << " " << BondCombinations[i][c];
        std::cout << "\n";
#endif
        auto resultEmbed = this->ComputeEmbeddingNumberCombo(container, this->BondCombinations[i]);

        if (this->Parameters.UseJonasFormat()) /// Jonas' format: g6 E(p_1) E(p_2) ... E(p_N,) where E(p_i) is the ith partition of the number of bonds based on the allowed bond lengths
        {
            if (i==0)
                fprintf(fpo, "%.*s", static_cast<int>(strlen(s)-1), s);
            fprintf(fpo, " %d", resultEmbed.first);
        }
        else /// default format
        {
            fprintf(fpo, "%.*s ", static_cast<int>(strlen(s)-1), s); /// write g6 string
            for (int d=0; d<this->MaxDegreeNeighbor; ++d) /// write bond combination i.e. parition of total number of bonds
                fprintf(fpo, "%d ", this->BondCombinations[i][d]);
            fprintf(fpo, "%d %d\n", resultEmbed.first, symmFactor);
            fflush(fpo);
        }
    }
    if (this->Parameters.UseJonasFormat())
    {
        fprintf(fpo, "\n");
        fflush(fpo);
    }
}

/// given a graph, compute embedding number for all bonds being bonds of the lattice (nearest-neighbors)
/// after calling ComputeEmbeddingNumberCombo, EmbedLists contains all of the embeddings
/// construct set of canonical graphs wrt cubic symmetry from EmbedLists and get counts of how many embeddings correspond to each canonical graph
/// @param container: GraphContainer object describing graph to be embedded
/// @param bondCombo: counts of types of bonds (NN, NNN, etc)
std::pair<std::vector<VertexEmbedList>, std::vector<int>> GraphEmbedder::ComputeCanonicalGraphsAndEmbeddingNumbers(GraphContainer container, const std::vector<int>& bondCombo)
{
    std::vector<VertexEmbedList> canonicalList; /// canonical list
    std::vector<VertexEmbedList> exampleList; /// one example from each canonical graph (return)
    std::vector<int> counts; /// counts of canonical lists to return

#ifdef DEBUG
    /// check that bond combo is valid
    if (this->MaxDegreeNeighbor!=bondCombo.size())
        throw std::invalid_argument("ComputeCanonicalGraphsAndEmbeddingNumbers requires bondCombo to be of size MaxDegreeNeighbor!\n");
    int result = 0;
    for (int i=0; i<bondCombo.size(); ++i)
        result += bondCombo[i];
    if (result!=container.GetL())
        throw std::invalid_argument("ComputeCanonicalGraphsAndEmbeddingNumbers requires bondCombo to be consistent with container!\n");
#endif

    auto countAndOneEmbedding = this->ComputeEmbeddingNumberCombo(container, bondCombo); /// do the embeddings

    /// now in EmbedLists we have the data that we need to get the list of canonical graphs and their counts
    for (auto it=this->EmbedLists.begin(); it!=this->EmbedLists.end(); ++it)
    {
        CubicLatticeCanonicalizor tempCanonicalizor(container, this->Lattice, *it); /// get canonical graph wrt octahedral group
        auto tempCanonical = tempCanonicalizor.GetCanonical();
        tempCanonical.SetNbrChoicesForFirstBond(it->GetNbrChoicesForFirstBond()); /// need to add this by hand to canonical!
        auto tempCanonicalOld = tempCanonicalizor.GetCanonicalOld();
#ifdef DEBUG
        if (tempCanonical!=tempCanonicalOld)
        {
            std::cout << "ERROR: OLD_AND_NEW_CANONICAL_DO_NOT_MATCH! " << container.GetG6String() << "\n";
            std::cout << "ORIGINAL:\n";
            std::cout << *it;
            std::cout << "NEW_CANONICAL:\n";
            std::cout << tempCanonical;
            std::cout << "OLD_CANONICAL:\n";
            std::cout << tempCanonicalOld;
        }
#endif
        auto tempIt = std::find(canonicalList.begin(), canonicalList.end(), tempCanonical);
        if (tempIt!=canonicalList.end()) /// if it is already in the list get index and increment counts
        {
            int tempIndex = std::distance(canonicalList.begin(), tempIt);
            counts[tempIndex]++;
        }
        else /// if it is not already in the list add canonical variant to canonicalList, original to exampleList, and add entry to counts
        {
            canonicalList.push_back(tempCanonical);
            exampleList.push_back(*it);
            counts.push_back(1);
        }
    }

#ifdef DEBUG
    int testTotal = 0;
#endif
    for (int i=0; i<canonicalList.size(); ++i) /// multiply counts by NbrChoicesForFirstBond
    {
#ifdef DEBUG
        std::cout << "DEBUG_ComputeEmbeddingNumbersAndCanonicalGraphs " << i << " \n";
        std::cout << canonicalList[i] << "\n";
#endif
        counts[i] *= canonicalList[i].GetNbrChoicesForFirstBond();
#ifdef DEBUG
        testTotal += counts[i];
#endif
    }

#ifdef DEBUG
    if (testTotal!=countAndOneEmbedding.first)
        throw std::logic_error("ComputeEmbeddingNumbersAndCanonicalGraphs requires the sum of the embedding numbers to equal the total!\n");
#endif

    return std::pair<std::vector<VertexEmbedList>, std::vector<int>>(exampleList, counts);

}

/// compute embedding number for a graph (rooted or unrooted); differs from old routine in that at every step, each list produces lists of higher order embedding new vertex connected with only ONE of previously embedded vertices
/// @param container: graph to embed
/// @param bondCombo: partition of N_{bonds} = \sum_i \tilde{N}_i, i=1,...,N_{deg}
/// @return pair consisting of the count (first) and an example embedding (second)
std::pair<int, VertexEmbedList> GraphEmbedder::ComputeEmbeddingNumberCombo(const GraphContainer& container, const std::vector<int> &bondCombo)
{

    auto lists = CreateInitialVertexEmbedLists(container, bondCombo);

    if (this->Parameters.EmbedCorrelator() && lists.size()==0) /// empty list for correlator i.e. rooted vertices
    {
#ifdef DEBUG
        std::cout << "COULD NOT EMBED ROOTED GRAPH! (ROOTED VERTICES CONNECTED BY A BOND BUT BOND NOT POSSIBLE)\n";
#endif
        return std::pair<int, VertexEmbedList>(0, VertexEmbedList(this->Parameters.GetMaxEmbeddingLength(), static_cast<MaxInteractionLength>(this->Parameters.GetCorrelatorLength()))); /// no embedding is possible
    }

    int vertexCount = 2; /// two vertices already embedded in initialization...
    while (vertexCount < container.GetN())
    {
        for (auto itLists=lists.begin(); itLists->GetSize()==vertexCount; ++itLists) /// ordered by number of number of vertices so can stop when we get to newly inserted VertexEmbedList objects
        {
            auto remainingVertices = this->GetRemainingVertices(*itLists); /// remaining vertices that have not been embedded
#ifdef DEBUG
            if (remainingVertices.size()!=(container.GetN()-vertexCount))
            {
                std::cerr << "Discrepancy with count and number of remaining vertices! " << remainingVertices.size() << " vs " << (container.GetN()-vertexCount) << "\n";
                std::exit(1);
            }
#endif
            auto dAllowedList = this->GetAllowedBondDegreesOfList(*itLists, bondCombo); /// list of allowed degrees based on bond combo and current VertexEmbedList object

            auto resultNextVertexToEmbed = this->DetermineNextVertexToEmbed(container, *itLists, remainingVertices); /// get vertex to embed next and vertex it is connected to among embedded vertices

            /// loop over allowed degrees of bonds
            for (int d=0; d<dAllowedList.size(); ++d)
            {
                for (int nn=1; nn<=this->GetNbrNeighbors(dAllowedList[d]); ++nn) /// try placing at dth-degree neighbor index nn of jth element of ith list
                {
#ifdef DEBUG
                    int tempIndexEnd = this->GetNeighbor(dAllowedList[d], resultNextVertexToEmbed.second.Index, nn);
                    std::vector<unsigned int> tempIndicesEnd(this->Lattice->GetDim(), 0);
                    this->Lattice->GetSiteCoordinates(tempIndexEnd, tempIndicesEnd);

                    std::vector<unsigned int> tempIndicesStart(this->Lattice->GetDim(), 0);
                    this->Lattice->GetSiteCoordinates(resultNextVertexToEmbed.second.Index, tempIndicesStart);
#endif
                    int nnIndex;
                    if (this->IsProposedNeighborSiteFree(*itLists, resultNextVertexToEmbed.second, dAllowedList[d], nn, nnIndex))
                    {
                        /// adding atleast one bond of degree d (will find out how many others)
                        std::vector<int> bondCountsToBeAdded(this->MaxDegreeNeighbor,0);
                        bondCountsToBeAdded[dAllowedList[d]] = 1;
                        if (this->IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(*itLists, container, bondCombo, resultNextVertexToEmbed.second, nnIndex, resultNextVertexToEmbed.first, bondCountsToBeAdded))
                        {
                            /// create new list and add to lists
                            VertexEmbedList temp(*itLists);
                            temp.AddVertexEmbed(VertexEmbed{resultNextVertexToEmbed.first,nnIndex}); /// add on new vertex to list

                            /// update bond counts
                            this->UpdateBondCounts(temp, bondCountsToBeAdded);

                            lists.insert(temp); /// automatically ordered

                        } /// if consistent

                    } /// if site free

                } /// for nn
            } /// for d

        } /// for itLists

        vertexCount++; /// increment vertex count

        /// check that all lists have appropriate number of elements and erase those that don't
        auto tempIt = lists.begin();
        while (tempIt!=lists.end())
        {
            if (tempIt->GetSize()!=vertexCount)
                tempIt = lists.erase(tempIt);
            else
                ++tempIt;
        }

        if (lists.size()==0) /// check if anything left?
            break;

    } /// while

    int result;
    if (this->Parameters.EmbedCorrelator())
    {
        result = lists.size(); /// no choices for placement of rooted vertices!
    }
    else
    {
        int temp = 0;
        for (auto it=lists.begin(); it!=lists.end(); ++it)
            temp += it->GetNbrChoicesForFirstBond();
        result = temp;
    }

#ifdef DEBUG
    std::cout << "FOUND " << lists.size() << " embeddings for graph!\n";
    std::cout << "THIS GIVES AN EMBEDDING NUMBER OF: " << result << "\n";

    int count = 0;
    for (auto itLists=lists.begin(); itLists!=lists.end(); ++itLists)
    {
        if (itLists->HasRepeatedSites())
            std::cerr << "ERROR! In ComputeEmbeddingNumberCombo embedding has repeated lattice sites!\n";
        if (itLists->HasRepeatedVertices())
            std::cerr << "ERROR! In ComputeEmbeddingNumberCombo embedding has repeated vertices!\n";
        std::cout << "Embedded graph number " << count;
        if (this->Parameters.EmbedCorrelator())
        {
            std::cout << " DISTANCE_OF_CORRELATOR " << itLists->GetCorrelatorDistanceAsIndex();
        }
        std::cout << "\n";
        for (auto itVertex=itLists->begin(); itVertex!=itLists->end(); ++itVertex)
        {
            std::vector<unsigned int> indices(this->Lattice->GetDim(), 0);
            this->Lattice->GetSiteCoordinates(itVertex->Index, indices);
            std::cout << "Vertex " << *itVertex;
            for (int k=0; k<this->Lattice->GetDim(); ++k)
                std::cout << " " << indices[k];
            std::cout << "\n";
        }
        count++;
    }
#endif
    this->EmbedLists = lists; /// save most recent list of embedded graphs
    if (lists.size()==0) /// check size
        return std::pair<int, VertexEmbedList>(0, VertexEmbedList(this->Parameters.GetMaxEmbeddingLength(), static_cast<MaxInteractionLength>(this->Parameters.GetCorrelatorLength())));
    else /// atleast one embedding to return
        return std::pair<int, VertexEmbedList>(result, VertexEmbedList(*lists.begin()));
}

/// compute embedding number for a graph (rooted or unrooted)
int GraphEmbedder::ComputeEmbeddingNumberComboOldWorking(const GraphContainer& container, const std::vector<int> &bondCombo)
{

    auto lists = CreateInitialVertexEmbedLists(container, bondCombo);

    if (this->Parameters.EmbedCorrelator() && lists.size()==0) /// empty list for correlator i.e. rooted vertices
    {
#ifdef DEBUG
        std::cout << "COULD NOT EMBED ROOTED GRAPH! (ROOTED VERTICES CONNECTED BY A BOND BUT BOND NOT POSSIBLE)\n";
#endif
        return 0; /// no embedding is possible
    }

    int vertexCount = 2; /// two vertices already embedded in initialization...
    while (vertexCount < container.GetN())
    {
        for (auto itLists=lists.begin(); itLists->GetSize()==vertexCount; ++itLists) /// ordered by number of number of vertices so can stop when we get to newly inserted VertexEmbedList objects
        {
            auto remainingVertices = this->GetRemainingVertices(*itLists); /// remaining vertices that have not been embedded
#ifdef DEBUG
            if (remainingVertices.size()!=(container.GetN()-vertexCount))
            {
                std::cerr << "Discrepancy with count and number of remaining vertices! " << remainingVertices.size() << " vs " << (container.GetN()-vertexCount) << "\n";
                std::exit(1);
            }
#endif
            auto dAllowedList = this->GetAllowedBondDegreesOfList(*itLists, bondCombo); /// list of allowed degrees based on bond combo and current VertexEmbedList object

            for (auto v=remainingVertices.begin(); v!=remainingVertices.end(); ++v) /// loop over vertices to be placed
            {
                for (auto itVertexEmbed=itLists->begin(); itVertexEmbed!=itLists->end(); ++itVertexEmbed)
                {
                    if (container.GetElementAdjacencyMatrix(*v,itVertexEmbed->Number))
                    {
                        /// loop over allowed degrees of bonds
                        for (int d=0; d<dAllowedList.size(); ++d)
                        {
                            for (int nn=1; nn<=this->GetNbrNeighbors(dAllowedList[d]); ++nn) /// try placing at dth-degree neighbor index nn of jth element of ith list
                            {
                                int nnIndex;
                                if (this->IsProposedNeighborSiteFree(*itLists, *itVertexEmbed, dAllowedList[d], nn, nnIndex))
                                {
                                    /// adding atleast one bond of degree d (will find out how many others)
                                    std::vector<int> bondCountsToBeAdded(this->MaxDegreeNeighbor,0);
                                    bondCountsToBeAdded[dAllowedList[d]] = 1;
                                    if (this->IsProposedSiteConsistentWithPreviousVerticesAndBondCounts(*itLists, container, bondCombo, *itVertexEmbed, nnIndex, *v, bondCountsToBeAdded))
                                    {
                                        /// create new list and add to lists
                                        VertexEmbedList temp(*itLists);
                                        temp.AddVertexEmbed(VertexEmbed{*v,nnIndex}); /// add on new vertex to list

                                        /// update bond counts
                                        this->UpdateBondCounts(temp, bondCountsToBeAdded);

                                        lists.insert(temp); /// insert (automatically ordered)

                                    } /// if consistent
                                } /// if site free

                            } /// for nn
                        } /// for d
                    } /// if edge exists

                } /// for j
            } /// for v

        } /// for i

        vertexCount++; /// increment vertex count

        /// check that all lists have appropriate number of elements and erase those that don't
        auto tempIt = lists.begin();
        while (tempIt!=lists.end())
        {
            if (tempIt->GetSize()!=vertexCount)
                tempIt = lists.erase(tempIt);
            else
                ++tempIt;
        }

        if (lists.size()==0) /// check if anything left?
            break;

    } /// while

    int result;
    if (this->Parameters.EmbedCorrelator())
    {
        result = lists.size(); /// no choices for placement of rooted vertices!
    }
    else
    {
        int temp = 0;
        for (auto it=lists.begin(); it!=lists.end(); ++it)
            temp += it->GetNbrChoicesForFirstBond();
        result = temp;
    }

#ifdef DEBUG
    std::cout << "FOUND " << lists.size() << " embeddings for graph!\n";
    std::cout << "THIS GIVES AN EMBEDDING NUMBER OF: " << result << "\n";

    int count = 0;
    for (auto itLists=lists.begin(); itLists!=lists.end(); ++itLists)
    {
        std::cout << "Embedded graph number " << count;
        if (this->Parameters.EmbedCorrelator())
        {
            std::cout << " DISTANCE_OF_CORRELATOR " << itLists->GetCorrelatorDistanceAsIndex();
        }
        std::cout << "\n";
        for (auto itVertex=itLists->begin(); itVertex!=itLists->end(); ++itVertex)
        {
            std::vector<unsigned int> indices(this->Lattice->GetDim(), 0);
            this->Lattice->GetSiteCoordinates(itVertex->Index, indices);
            std::cout << "Vertex " << *itVertex;
            for (int k=0; k<this->Lattice->GetDim(); ++k)
                std::cout << " " << indices[k];
            std::cout << "\n";
        }
        count++;
    }
#endif

   return result;
}

/// update the counts of each type of bonds for a list of embedded vertices (called each time we embed a new vertex)
/// @param list: list of embedded vertices
/// @param countsToBeAdded: vector giving the counts for each type of bond (NN, NNN, etc.)
void GraphEmbedder::UpdateBondCounts(VertexEmbedList &List, const std::vector<int>& countsToBeAdded)
{
    if (countsToBeAdded.size()!=this->MaxDegreeNeighbor)
        throw std::invalid_argument("UpdateBondCounts requires countsToBeAdded to be of size MaxDegreeNeighbor!\n");
    for (int b=0; b<countsToBeAdded.size(); ++b)
        for (int k=0; k<countsToBeAdded[b]; ++k)
            List.IncrementBondCount(b);
}

void GraphEmbedder::PrintVertexEmbedList(const VertexEmbedList& list)
{
    std::vector<unsigned int> tempIndices(this->Lattice->GetDim(), 0);
    std::vector<unsigned int> indicesOrigin(this->Lattice->GetDim(), this->Lattice->GetN()/2);

    std::cout << "********VertexEmbedList********\n";
    for (auto it=list.begin(); it!=list.end(); ++it)
    {
        std::cout << "Vertex " << it->Number << " at site: (shifted)";
        this->Lattice->GetSiteCoordinates(it->Index, tempIndices);
        for (auto i=0; i<tempIndices.size(); ++i)
            std::cout << " " << int(tempIndices[i])-int(indicesOrigin[i]);
        std::cout << " (unshifted)";
        for (auto i=0; i<tempIndices.size(); ++i)
            std::cout << " " << tempIndices[i];
        std::cout << "\n";
    }
}


/// find first nonzero element in the adjacency matrix (search from the "back") and then place vertices at first NN
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

void GraphEmbedder::GetCombinationsOfBondsFixedManhattanDistance(int nbrBonds, int manhattanDistance)
{
    if (manhattanDistance<=0)
        throw std::invalid_argument("GetCombinationsOfBondsFixedManhattanDistance requires a positive value for argument manhattanDistance!\n");

    std::vector<int> bondCounts(this->MaxDegreeNeighbor,-1);
    std::vector<int> List(nbrBonds+1);
    for (int i = 0; i<List.size(); ++i)
        List[i] = i;

    std::vector<int> bondWeightsManhattan(this->MaxDegreeNeighbor,-1);
    for (int i=0; i<bondWeightsManhattan.size(); ++i)
        bondWeightsManhattan[i] = this->Lattice->GetManhattanDistance(i+1);

    auto correctManhattan = [bondWeightsManhattan,manhattanDistance](const std::vector<int>& partition) /// lambda function
    {
        if (bondWeightsManhattan.size()!=partition.size())
            throw std::invalid_argument("Bond weights in capture should be the same size as partition!\n");
        return (std::inner_product(bondWeightsManhattan.begin(),bondWeightsManhattan.end(),partition.begin(),0.0)==manhattanDistance);
    };

    this->BondCombinations.clear(); /// clear list of bonds

    this->GenerateCombinations(List, bondCounts, 0, correctManhattan); /// generate combinations

}

/// get all combinations of counts for bond types {N_i, i=0(NN),1(NNN),etc. } where \sum_i N_i = nbrBonds
/// nbrBonds: total number of bonds for a particular graph
void GraphEmbedder::GetCombinationsOfBondsFixedNumberOfBonds(int nbrBonds)
{
    std::vector<int> bondCounts(this->MaxDegreeNeighbor,-1);
    std::vector<int> List(nbrBonds+1);
    for (int i = 0; i<List.size(); ++i)
        List[i] = i;

    auto partitionEqualsBonds = [nbrBonds](const std::vector<int>& partition)
    {
        return (std::accumulate(partition.begin(), partition.end(), 0)==nbrBonds);
    };

    this->BondCombinations.clear(); /// clear list of bonds

    this->GenerateCombinations(List, bondCounts, 0, partitionEqualsBonds); /// generate combinations

}

/// get all combinations of vertices using depth first search
void GraphEmbedder::GenerateCombinations(const std::vector<int>& arr, std::vector<int>& data, int index, std::function<bool(const std::vector<int>&)> isValid)
{

    if (index == data.size()) /// created a combo?
    {
        if (isValid(data))
        {
            this->BondCombinations.push_back(data);
        }
        return;
    }

    for (int i=0; i<arr.size(); ++i)
    {
        data[index] = arr[i];
        this->GenerateCombinations(arr, data, index+1, isValid);
    }

}

/// calls the appropriate routine for the initial embedded vertices
/// @param container: graph to be embedded
/// @param bondCombo: partition of N_{bonds} to embed
std::set<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedLists(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    if (this->Parameters.EmbedCorrelator())
        return CreateInitialVertexEmbedListsRooted(container, bondCombo);
    else
        return CreateInitialVertexEmbedListsNonRooted(container, bondCombo);
}

/// creates initial list of embedded vertices for rooted graphs
/// NOTE: embedded vertices ALWAYS taken to be labeled 1 and 2 (0 and 1 in NAUTY convention)! Convention specified in GraphGeneratorNauty
std::set<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRooted(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    std::vector<int> rootedVertices{1,2}; /// CONVENTION!
    return this->CreateInitialVertexEmbedListsRootedFixed(container, bondCombo, rootedVertices);
}

/// create initial list for the color given a pair of fixed vertices and a set of bond specifications
/// @param container: graph to be embedded
/// @param bondCombo: set of bond specifications
/// @param rootedVertices: labels of rooted vertices
std::set<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsRootedFixed(const GraphContainer& container, const std::vector<int> &bondCombo, const std::vector<int>& rootedVertices)
{
    std::set<VertexEmbedList> result;

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
        result.insert(tempList);
    }
#ifdef DEBUG
    std::cout << "INITIAL_LIST:\n";
    for (auto it=result.begin(); it!=result.end(); ++it)
        std::cout << " " << *it << "\n";
#endif
    return result;
}

/// creates the initial list for non-rooted graphs
/// looks for the first non-zero off-diagonal element of the adjacency matrix and then embeds that bond on the lattice in the canonically defined "first direction" for all allowed bond lengths
/// @param container: graph to be embedded
/// @param bondCombo: set of bond counts
std::set<VertexEmbedList> GraphEmbedder::CreateInitialVertexEmbedListsNonRooted(const GraphContainer& container, const std::vector<int> &bondCombo)
{
    std::set<VertexEmbedList> result;

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
                    result.insert(tempList);
                }
            }
            break;
        }
    }
    return result;
}

/// from g6 string in command line arguments, compute the canonical embeddings (NN only) on the cubic lattice and their counts
/// container: graph container consistent with g6 string given to parameters (this should come from CanonicalGraphManager!!!!!)
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> GraphEmbedder::ComputeCanonicalEmbeddingsAndCountsNN(const GraphContainer &container)
{
    if (this->Parameters.GetG6().empty())
        throw std::invalid_argument("GraphEmbedder::ComputeCanonicalGraphsAndCountsNN requires the user to input a g6 string!\n");

    if (this->ResolveLatticeType(this->Parameters.GetLatticeType())!=GraphEmbedder::LatticeType::Cubic)
        throw std::invalid_argument("GraphEmbedder::ComputeCanonicalGraphsAndCountsNN requires the cubic lattice!\n");

    if (container.GetG6String()!=this->Parameters.GetG6()) /// make sure things are consistent
        throw std::invalid_argument("GraphEmbedder::ComputeCanonicalGraphsAndCountsNN requires container to be consistent with g6 string!\n");

    if (this->Parameters.EmbedCorrelator() && container.GetNbrRooted()!=2) /// check consistency
        throw std::invalid_argument("GraphEmbedder::ComputeCanonicalGraphsAndCountsNN requires container to be two-rooted if embedding a correlators!\n");

    if (!this->Parameters.EmbedCorrelator() && container.GetNbrRooted()!=0) /// check consistency
        throw std::invalid_argument("GraphEmbedder::ComputeCanonicalGraphsAndCountsNN requires container to be unrooted if not embedding a correlator!\n");

    /// all bonds are NN i.e. true embedding on cubic lattice
    std::vector<int> bondCombo(this->MaxDegreeNeighbor, 0);
    bondCombo[0] = container.GetL();

    /// get canonical graphs and counts
    auto canonicalListAndCounts = this->ComputeCanonicalGraphsAndEmbeddingNumbers(container, bondCombo);
#ifdef DEBUG
    std::cout << "SIZES_DEBUG: " << canonicalListAndCounts.first.size() << " " << canonicalListAndCounts.second.size() << "\n";
    for (int i=0; i<canonicalListAndCounts.first.size(); ++i)
    {
        std::cout << "embedding " << i << " with counts " << canonicalListAndCounts.second[i] << "\n";
        std::cout << canonicalListAndCounts.first[i] << "\n";
    }
#endif

    /// combine container with graphs and counts
    return std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>>(container, canonicalListAndCounts.first, canonicalListAndCounts.second);

}

/// do not use with rooted graphs! each simple connected graph generates some number of rooted graphs and they are written to file
/// quit if we ask for rooted graphs as we are unsure how to streamline the procedure from the files????
std::pair<GraphContainer, VertexEmbedList> GraphEmbedder::ContainerAndSampleCubicEmbeddingFromG6()
{
    if (this->Parameters.EmbedCorrelator())
        throw std::invalid_argument("GraphEmbedder::ContainerFromG6 not embedding rooted graphs!");

    if (this->Parameters.GetG6().empty())
        throw std::invalid_argument("GraphEmbedder::ContainerFromG6 requires the user to input a g6 string!\n");

    if (this->ResolveLatticeType(this->Parameters.GetLatticeType())!=GraphEmbedder::LatticeType::Cubic)
        throw std::invalid_argument("GraphEmbedder::ContainerFromG6 requires the cubic lattice!\n");

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc"); /// allocate graph

    this->DenseNautyFromString(this->G6String, g);

    GraphContainer container(this->N, this->MWords, g); /// container from densenauty

    this->GetCombinationsOfBondsFixedNumberOfBonds(container.GetL()); /// generate combos of bonds for fixed bond number

    auto resultEmbed = this->ComputeEmbeddingNumberCombo(container, this->BondCombinations[0]);

    if (resultEmbed.first==0)
        std::cerr << "ERROR in ContainerAndSampleCubicEmbeddingFromG6(): graph has no embeddings on the cubic lattice!\n";

    DYNFREE(g, g_sz); /// free graph

    return std::pair<GraphContainer, VertexEmbedList>(container, resultEmbed.second);
}

/// get the number of vertices from a g6 string
int GraphEmbedder::GetGraphSizeFromString(const std::string& g6String)
{
    char *tempg6 = new char[g6String.length()+1];
    std::strcpy(tempg6, g6String.c_str());
    int result = graphsize(tempg6);
    delete[] tempg6;
    return result;
}

/// get the densenauty form of a graph from the g6 string
void GraphEmbedder::DenseNautyFromString(const std::string& g6String, graph *g)
{
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, this->MWords); /// g6 string to densenauty

    delete[] tempg6;
}

/// get the symmetry factor from a given unrooted graph
/// g: pointer to allocated graph of size N
int GraphEmbedder::GetSymmFactor(graph *g)
{
    DYNALLSTAT(int,lab,lab_sz); /// declare labels
    DYNALLSTAT(int,ptn,ptn_sz); /// declare pattern
    DYNALLSTAT(int,orbits,orbits_sz); /// declare orbits
    statsblk status; /// status

    DYNALLOC2(int, lab, lab_sz, this->N, this->MWords, "malloc"); /// allocate labels
    DYNALLOC2(int, ptn, ptn_sz, this->N, this->MWords, "malloc"); /// allocate pattern
    DYNALLOC2(int, orbits, orbits_sz, this->N, this->MWords, "malloc"); /// allocate orbits

    this->CallDenseNauty(g, lab, ptn, orbits, status); /// call NAUTY

    /// write out to temporary file and read back in
    std::string tempFilename = "temp_symm"; /// HARD-CODED
    FILE *fpo = fopen(tempFilename.c_str(), "w"); /// open file for writing
    if (fpo!=NULL)
        std::cout << "GraphEmbedder::GetSymmFactor: Successfully opened " << tempFilename << " for writing!\n";
    else
        throw std::invalid_argument("GraphEmbedder::GetSymmFactor: Error opening "+tempFilename);

    writegroupsize(fpo,status.grpsize1,status.grpsize2); /// write symmetry factor to file
    fprintf(fpo, "\n");

    fclose(fpo); /// close file

    FILE *fpi = fopen(tempFilename.c_str(), "r"); /// open file for reading
    if (fpi!=NULL)
        std::cout << "GraphEmbedder::GetSymmFactor: Successfully opened " << tempFilename << " for reading!\n";
    else
        throw std::invalid_argument("GraphEmbedder::GetSymmFactor: Error opening "+tempFilename);

    int symmFactor = -1;
    fscanf(fpi, "%d\n", &symmFactor); /// read symmetry factor
    if (symmFactor==-1)
        throw std::invalid_argument("GetSymmFactor symmFactor has wrong value!\n");

    fclose(fpi); /// close file

    DYNFREE(lab, lab_sz); /// free labels
    DYNFREE(ptn, ptn_sz); /// free pattern
    DYNFREE(orbits, orbits_sz); /// free orbits

    return symmFactor;
}

/// we require one option if the other is not provided by the user (option either default or not specified)
/// @param vm: variable map
/// @param required: name of the required option
/// @param other: the option which determines whether the first is actually required
void GraphEmbedderParametersNauty::RequiredOptionWhenOtherOptionMissing(const po::variables_map& vm, const char* required, const char* other)
{
    if(vm.count(other)==0 || vm[other].defaulted())
        if (vm.count(required) == 0 || vm[required].defaulted())
            throw std::logic_error(std::string("Option ")+required+" required when "+other+" missing!\n");
}

/// process user's command-line arguments
bool GraphEmbedderParametersNauty::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                ("order,n", po::value<unsigned int>(&this->N)->required(), "Vertex order")
                ("input,i", po::value<std::string>(&this->InputFilename), "Input filename")
                ("output,o", po::value<std::string>(&this->OutputFilename), "Output filename")
                ("g6,g", po::value<std::string>(&this->G6)->default_value(""), "Embed  single graph from g6 string")
                (",l", po::value<std::string>(&this->LatticeType)->default_value("Cubic"), "Lattice type: Cubic, Square, or Triangular")
                (",d", po::value<MaxInteractionLength>(&this->MaxEmbeddingLength)->default_value(MaxInteractionLength::NearestNeighbor), "MaxEmbeddingLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                (",c", po::bool_switch(&this->Correlator)->default_value(false), "Embed correlators (two-point functions)?")
                (",j", po::bool_switch(&this->JonasFormat)->default_value(false), "Output using Jonas' format?")
                (",m", po::value<MaxInteractionLength>(&this->CorrelatorLength)->default_value(MaxInteractionLength::NearestNeighbor), "CorrelatorLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        this->RequiredOptionWhenOtherOptionMissing(vm, "input", "g6");
        this->RequiredOptionWhenOtherOptionMissing(vm, "output", "g6");

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
