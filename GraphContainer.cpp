#include "GraphContainer.h"

#include <iostream>

/// constructor which accepts size and number of words (need to call SetGraphFromDenseNauty before container can be used!)
GraphContainer::GraphContainer(int n, int m) :
    N(n),
    MWords(m),
    L(0),
    NTimesNMinusOneDiv2(n*(n-1)/2),
    VertexOrder(n,0),
    M(n, std::vector<bool>(n,false)),
    RowM(n*(n-1)/2),
    ColM(n*(n-1)/2)
{
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
void GraphContainer::GetDenseNautyFromGraph(graph *g)
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
GraphContainer::GraphContainer(int n, int m, graph *g) :
    N(n),
    MWords(m),
    L(0),
    NTimesNMinusOneDiv2(n*(n-1)/2),
    VertexOrder(n,0),
    M(n, std::vector<bool>(n,false)),
    RowM(n*(n-1)/2),
    ColM(n*(n-1)/2)
{

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
/// @param labOld: N element array with initial labels (should match already set up adjacency matrix) (starts at zero!)
/// @param labNew: N element array with final labels (starts at zero!)
/// @param v1: label of first rooted vertex (starts at zero!) (debugging: make sure map takes v1 to v1)
/// @param v2: label of second rooted vertex (starts at zero!) (debugging: make sure map takes v2 to v2)
void GraphContainer::ColoredCanonicalRelabeling(int *labOld, int *labNew, int v1, int v2, bool verbose)
{
    std::vector<int> alpha(this->N); /// map of isomorphism

    for (int i=0; i<this->N; ++i)
    {
        alpha[labOld[i]] = labNew[i];
    }

    if (!(alpha[v1]==v1 && alpha[v2]==v2))
        throw std::invalid_argument("ColoredCanonicalRelabeling v1 and v2 both must map to themselves!\n");

    if (verbose)
    {
        for (int i=0; i<this->N; ++i)
            std::cout << "alpha maps vertex " << i << " to vertex " << alpha[i] << "\n";
    }

    std::vector<std::vector<bool>> newM(N, std::vector<bool>(this->N, false));
    for (int i=0; i<this->N; ++i)
    {
        for (int j=i+1; j<this->N; ++j)
        {
            newM[i][j] = this->GetElementAdjacencyMatrix(alpha[i]+1, alpha[j]+1); /// RowM and ColM start at 1 and alpha at 0!!!!
            newM[j][i] = newM[i][j]; /// symmetrize
            //std::cout << "DEBUG_NEWM: " << i << " " << this->ColM[i] << " " << newM[this->RowM[i]-1][this->ColM[i]-1] << "\n";
        }
    }

    this->M = newM; /// copy adjacency matrix
    /// set the vertex orders
    for (unsigned int v=1; v<=this->N; ++v)
        this->SetVertexOrder(v, this->ComputeVertexOrder(v));

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

/// compare two containers
/// compare number of vertices, number of bonds and finally adjacency matrices
bool operator==(const GraphContainer& lhs, const GraphContainer& rhs)
{
    if (lhs.GetN()!=rhs.GetN()) /// compare graph order
        return false;
    if (lhs.GetL()!=rhs.GetL()) /// compare number of bonds
        return false;
    for (int i=1; i<lhs.GetN(); ++i) /// compare elements of adjacency matrix
        for (int j=0; j<i; ++j)
            if (lhs.GetElementAdjacencyMatrix(i+1,j+1)!=rhs.GetElementAdjacencyMatrix(i+1,j+1))
                return false;
    return true;
}

/// inequality uses the result of comparison
bool operator!=(const GraphContainer& lhs, const GraphContainer& rhs)
{
    return !(lhs==rhs);
}
