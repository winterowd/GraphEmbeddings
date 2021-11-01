#include "GraphContainer.h"

#include <iostream>

/// constructor which accepts size and number of words (need to call SetGraphFromDenseNauty before container can be used!)
GraphContainer::GraphContainer(int n, int m, bool storeRooted, int nbrRooted) :
    N(n),
    MWords(m),
    L(0),
    StoreRooted(storeRooted),
    NbrRooted(nbrRooted),
    RootedVertices(nbrRooted,-1),
    NTimesNMinusOneDiv2(n*(n-1)/2),
    VertexOrder(n,0),
    M(n, std::vector<bool>(n,false)),
    RowM(n*(n-1)/2),
    ColM(n*(n-1)/2)
{   
    if (storeRooted && (nbrRooted<0 || nbrRooted>2))
        throw std::invalid_argument("GraphContainer constructor needs 0 < nbrRooted <= 2 if storeRooted is true!\n");

    this->SetRowAndColMapping();
}

/// mapping from linear index traversing upper triangular part of adjacency matrix to rows and columns
void GraphContainer::SetRowAndColMapping()
{
    /// fill look up tables for upper triangular elements of adjacency matrix
    int b=-1;
    for (int i=1; i<this->N; ++i)
    {
        for (int j=0; j<i; ++j)
        {
            b++;
            this->RowM[b] = j+1;
            this->ColM[b] = i+1;
        }
    }
}

/// convert graph to dense NAUTY structure
/// g: pointer to dense nauty data structure (memory assumed to be allocated elsewhere!)
void GraphContainer::GetDenseNautyFromGraph(graph *g) const
{
    if (g==NULL)
        throw std::invalid_argument("SetDenseNautyFromGraph expects g to be point to memory that is already allocated!\n");

    EMPTYGRAPH(g, this->N, this->MWords); /// clear graph

    for (int j=1; j<this->N; ++j)
    {
        for (int i=0; i<j; ++i)
        {
            if (this->GetElementAdjacencyMatrix(i+1,j+1))
                ADDONEEDGE(g, i, j, this->MWords);
        }
    }
}

/// set adjacency matrix and all other private variables from dense nauty data structure
/// NOTE: assumes N and MWords set!
/// g: pointer to dense nauty data structure
void GraphContainer::SetGraphFromDenseNauty(graph *g)
{
    set *gj;
    auto count = 0;
    for (int j=1; j<this->N; ++j)
    {
        gj = GRAPHROW(g,j,this->MWords);
        for (int i=0; i<j; ++i)
        {
            if (ISELEMENT(gj,i))
            {
                this->SetElementAdjacenyMatrix(i+1, j+1);
                this->SetElementAdjacenyMatrix(j+1, i+1);
                count++;
            }
            else
            {
                this->UnsetElementAdjacenyMatrix(i+1, j+1);
                this->UnsetElementAdjacenyMatrix(j+1, i+1);
            }
        }
    }

    this->L = count; /// set number of bonds

    /// set the vertex orders
    for (unsigned int v=1; v<=this->N; ++v)
        this->SetVertexOrder(v, this->ComputeVertexOrder(v));

#ifdef DEBUG
    this->PrintM();
    this->PrintVertexOrders();
#endif

}

/// constructor from dense nauty data structure
GraphContainer::GraphContainer(int n, int m, graph *g, bool storeRooted, int nbrRooted) :
    N(n),
    MWords(m),
    L(0),
    StoreRooted(storeRooted),
    NbrRooted(nbrRooted),
    RootedVertices(nbrRooted,-1),
    NTimesNMinusOneDiv2(n*(n-1)/2),
    VertexOrder(n,0),
    M(n, std::vector<bool>(n,false)),
    RowM(n*(n-1)/2),
    ColM(n*(n-1)/2)
{
    if (storeRooted && (nbrRooted<0 || nbrRooted>2))
        throw std::invalid_argument("GraphContainer constructor needs 0 < nbrRooted <= 2 if storeRooted is true!\n");

    this->SetRowAndColMapping();

    this->SetGraphFromDenseNauty(g);

}

/// relabel vertices of the graph
/// newLabels: N element vector with the new labels i.e. newLabels_i = j means vertex i mapped to vertex j (starts at zero!)
void GraphContainer::RelabelVertices(const std::vector<int>& newLabels)
{
    if (newLabels.size()!=this->N)
        throw std::invalid_argument("RelabelVertices requires newLabels to be of size N!\n");

    std::vector<std::vector<bool>> newM(N, std::vector<bool>(this->N, false));
    for (int i=0; i<this->N; ++i)
    {
        for (int j=i+1; j<this->N; ++j)
        {
            newM[i][j] = this->GetElementAdjacencyMatrix(newLabels[i]+1, newLabels[j]+1); /// newLabels assumes vertex labels starting at zero!
            newM[j][i] = newM[i][j]; /// symmetrize
        }
    }

    this->M = newM; /// copy adjacency matrix
    /// set the vertex orders
    for (unsigned int v=1; v<=this->N; ++v)
        this->SetVertexOrder(v, this->ComputeVertexOrder(v));

#ifdef DEBUG
    std::cout << "AFTER_RELABELING:\n";
    this->PrintM();
    this->PrintVertexOrders();
#endif
}

/// for two-rooted (colored) graphs, after calling densenauty to get the canonical form
/// @param labCanon: N element array giving the mapping to the canonical graph i.e. vertex labCanon[i] is relabeled as vertex i
/// @param v1: label of first rooted vertex (starts at zero!) (debugging: make sure map takes v1 to v1)
/// @param v2: label of second rooted vertex (starts at zero!) (debugging: make sure map takes v2 to v2)
void GraphContainer::ColoredCanonicalRelabeling(int *labCanon, int v1, int v2, bool verbose)
{

    if (labCanon[0]!=v1)
        throw std::invalid_argument("ColoredCanonicalRelabeling labCanon[0] must equal v1!\n");

    if (v2!=-1 && labCanon[1]!=v2)
        throw std::invalid_argument("ColoredCanonicalRelabeling labCanon[1] must equal v2!\n");

    if (verbose)
    {
        for (int i=0; i<this->N; ++i)
            std::cout << "labCanon maps vertex " << labCanon[i] << " to vertex " << i << "\n";
    }

    std::vector<std::vector<bool>> newM(N, std::vector<bool>(this->N, false));
    for (int i=0; i<this->N; ++i)
    {
        for (int j=i+1; j<this->N; ++j)
        {
            newM[i][j] = this->GetElementAdjacencyMatrix(labCanon[i]+1,labCanon[j]+1);  /// RowM and ColM start at 1 and alpha at 0!!!!
            newM[j][i] = newM[i][j]; /// symmetrize
        }
    }

    this->M = newM; /// copy adjacency matrix
    /// set the vertex orders
    for (unsigned int v=1; v<=this->N; ++v)
        this->SetVertexOrder(v, this->ComputeVertexOrder(v));

    if (this->StoreRooted)
    {
        if (v1!=this->RootedVertices[0])
            throw std::invalid_argument("ColoredCanonicalRelabeling requires RootedVertices[0] to equal v1!\n");
        this->RootedVertices[0] = 0; /// first rooted vertex is always labeled zero
        if (v2!=-1 && v2!=this->RootedVertices[1])
            throw std::invalid_argument("ColoredCanonicalRelabeling requires RootedVertices[1] to equal v2!\n");
        if (v2!=-1) /// check if second vertex is set
            this->RootedVertices[1] = 1; /// second rooted vertex is always labeled zero
    }

#ifdef DEBUG
    std::cout << "AFTER_COLORED_CANONICAL_RELABELING:\n";
    this->PrintM();
    this->PrintVertexOrders();
#endif

}

/// print out adjacency matrix
void GraphContainer::PrintM() const
{
    for (int i=0; i<this->NTimesNMinusOneDiv2; ++i)
        std::cout << "M: " << i << " " << RowM[i] << " " << ColM[i] << " " << this->GetElementAdjacencyMatrix(RowM[i], ColM[i]) << "\n";
}

/// accessor
int GraphContainer::GetRowM(int index) const
{
    if (index<0 || index>=this->NTimesNMinusOneDiv2)
        throw std::invalid_argument("GetRowM must have 0 <= index < NTimesNMinusOneDiv2!\n");
    return this->RowM[index];
}


/// accessor
int GraphContainer::GetColM(int index) const
{
    if (index<0 || index>=this->NTimesNMinusOneDiv2)
        throw std::invalid_argument("GetColM must have 0 <= index < NTimesNMinusOneDiv2!\n");
    return this->ColM[index];
}

/// print out vertex orders
void GraphContainer::PrintVertexOrders() const
{
    for (int i=1; i<=this->N; ++i)
        std::cout << "VERTEX_ORDER: " << i << " " << this->GetVertexOrder(i) << "\n";
}

/// check if vertex orders are consistent with current state of adjacency matrix
bool GraphContainer::VertexOrdersConsistentWithAdjacencyMatrix()
{
    auto result = true;
    for (int v=1; v<=this->N; ++v)
    {
        if (this->GetVertexOrder(v) != this->ComputeVertexOrder(v))
        {
            std::cerr << "ERROR! In VertexOrdersConsistentWithAdjacencyMatrix() sum over row for " << v << " does not agree with VertexOrder!" << "\n";
            result = false;
        }
    }
    return result;
}

/// check that adjacency matrix has zeros on the diagonal and is symmetric
bool GraphContainer::AdjacencyMatrixOK()
{
    auto result = true;
    for (int v1=1; v1<=this->N; ++v1)
    {
        if (this->GetElementAdjacencyMatrix(v1,v1))
        {
            std::cerr << "ERROR! In AdjacencyMatrixOK() diagonal element " << v1 << " is nonzero!" << "\n";
            result = false;
        }
        for (int v2=v1+1; v2<=this->N; ++v2)
        {
            if (this->GetElementAdjacencyMatrix(v1,v2) != this->GetElementAdjacencyMatrix(v2,v1))
            {
                std::cerr << "ERROR! In AdjacencyMatrixOK() element " << v1 << " " << v2 << " is not equal to element " << v2 << " " << v1 << "\n";
                result = false;
            }
        }
    }
    return result;
}

/// compute the vertex order from a row of the adjacency matrix
int GraphContainer::ComputeVertexOrder(unsigned int v)
{
    int result = 0;
    for (int i=0; i<this->M[v-1].size(); ++i) /// sum along row of adjacency matrix
        result += this->M[v-1][i];
    return result;
}

/// accessor for adjacency marix
bool GraphContainer::GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("GetElementAdjacencyMatrix 1 <= v1, v2 <= N!");
    return this->M[v1-1][v2-1];
}

/// accessor for adjacency marix
bool GraphContainer::GetElementAdjacencyMatrix(unsigned index) const
{
    if (index >= this->NTimesNMinusOneDiv2)
        throw std::invalid_argument("GetElementAdjacencyMatrix 0 <= index <= N(N-1)/2!");
    return this->GetElementAdjacencyMatrix(this->RowM[index], this->ColM[index]);
}

/// set element of adjacency matrix to 1
void GraphContainer::SetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("SetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    this->M[v1-1][v2-1] = true;
}

/// set element of adjacency matrix to 0
void GraphContainer::UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("UnsetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    this->M[v1-1][v2-1] = false;
}

/// accessor for vertex order
int GraphContainer::GetVertexOrder(unsigned int v) const
{
    if (v > this->N || v < 1)
        throw std::invalid_argument("GetVertexOrder requires N >= v >= 1!");
    return this->VertexOrder[v-1];
}

/// set the vertex order
void GraphContainer::SetVertexOrder(unsigned int v, int order)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("SetVertexOrder requires N >= v >= 1!");
    if (order > this->N || order < 0)
        throw std::invalid_argument("SetVertexOrder requires N >= order >= 1!");
    this->VertexOrder[v-1] = order;
}

/// set one of the rooted vertices to have a label
/// @param index: index of rooted vertices (color)
/// @param label: label of given rooted vertex
void GraphContainer::SetRootedVertex(int index, int label)
{
    if (!this->StoreRooted)
        std::cerr << "ERROR: SetRootedVertex called when StoreRooted flag is false!\n";
    if (index<0 || index>=this->NbrRooted)
        throw std::invalid_argument("SetRootedVertex requires 0 <= index < NbrRooted!\n");
    if (label<0 || label>=this->N)
        throw std::invalid_argument("SetRootedVertex requires 0 <= label < N!\n");
    this->RootedVertices[index] = label;
}

/// accessor to obtain label of rooted vertex of index
int GraphContainer::GetRootedVertex(int index) const
{
    if (!this->StoreRooted)
        std::cerr << "ERROR: GetRootedVertex called when StoreRooted flag is false!\n";
    if (index<0 || index>=this->NbrRooted)
        throw std::invalid_argument("GetRootedVertex requires 0 <= index < NbrRooted!\n");
    return this->RootedVertices[index];
}

/// accessor for number of rooted vertices
int GraphContainer::GetNbrRooted() const
{
    if (!this->StoreRooted)
        std::cerr << "ERROR: GetNbrRooted called when StoreRooted flag is false!\n";
    return this->NbrRooted;
}

/// construct bit string from adjacency matrix and return this as an integer
int GraphContainer::ComputeCurrentKey() const
{
    int result = 0;
    int count = 0;
    for (int i=this->NTimesNMinusOneDiv2-1; i>=0; --i)
    {
       result += this->GetElementAdjacencyMatrix(this->RowM[i],this->ColM[i]) << count;
       count++;
    }
    return result;
}

/// compare two containers
/// compare number of vertices, number of bonds and finally adjacency matrices
bool operator==(const GraphContainer& lhs, const GraphContainer& rhs)
{
    if (lhs.StoringRooted()!=rhs.StoringRooted())
        throw std::invalid_argument("ERROR: Comparing two containers and find that one is storing rooted vertices and the other is not!\n");
    //std::cout << "HI_INSIDE!\n";
    if (lhs.StoringRooted())
    {
        //std::cout << "HI_STORING_ROOTED!\n";
        if (lhs.GetNbrRooted()!=rhs.GetNbrRooted())
        {
            //std::cout << "DEBUG_EQUALITY_CONTAINER: DIFFERENT_NUMBER_ROOTED!\n";
            return false;
        }
        for (int i=0; i<lhs.GetNbrRooted(); ++i)
            if (lhs.GetRootedVertex(i)!=rhs.GetRootedVertex(i))
            {
                //std::cout << "DEBUG_EQUALITY_CONTAINER: DIFFERENT_ROOTED_VERTICES!\n";
                return false;
            }
        //std::cout << "HI_STORING_ROOTED_END!\n";
    }

    if (lhs.GetN()!=rhs.GetN()) /// compare graph order
    {
        //std::cout << "DEBUG_EQUALITY_CONTAINER: DIFFERENT_ORDER!\n";
        return false;
    }
    if (lhs.GetL()!=rhs.GetL()) /// compare number of bonds
    {
        //std::cout << "DEBUG_EQUALITY_CONTAINER: DIFFERENT_NBR_BONDS!\n";
        return false;
    }
    for (int i=1; i<lhs.GetN(); ++i) /// compare elements of adjacency matrix
        for (int j=0; j<i; ++j)
            if (lhs.GetElementAdjacencyMatrix(i+1,j+1)!=rhs.GetElementAdjacencyMatrix(i+1,j+1))
            {
                //std::cout << "DEBUG_EQUALITY_CONTAINER: DIFFERENT_ADJACENCY_MATRICES!\n";
                return false;
            }
    return true;
}

/// inequality uses the result of comparison
bool operator!=(const GraphContainer& lhs, const GraphContainer& rhs)
{
    return !(lhs==rhs);
}

bool operator<(const GraphContainer& lhs, const GraphContainer& rhs)
{
    return lhs.ComputeCurrentKey() < rhs.ComputeCurrentKey();
}
