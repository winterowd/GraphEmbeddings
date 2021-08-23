#include "GraphContainer.h"

#include <iostream>

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

/// NOTE: assumes N and MWords set!
void GraphContainer::SetGraphFromNauty(graph *g)
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

    this->SetGraphFromNauty(g);

}

void GraphContainer::PrintM() const
{
    for (int i=0; i<this->NTimesNMinusOneDiv2; ++i)
        std::cout << "M: " << i << " " << RowM[i] << " " << ColM[i] << " " << this->GetElementAdjacencyMatrix(RowM[i], ColM[i]) << "\n";
}

void GraphContainer::PrintVertexOrders() const
{
    for (int i=1; i<=this->N; ++i)
        std::cout << "VERTEX_ORDER: " << i << " " << this->GetVertexOrder(i) << "\n";
}

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

int GraphContainer::ComputeVertexOrder(unsigned int v)
{
    int result = 0;
    for (int i=0; i<this->M[v-1].size(); ++i) /// sum along row of adjacency matrix
        result += this->M[v-1][i];
    return result;
}

bool GraphContainer::GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("GetElementAdjacencyMatrix 1 <= v1, v2 <= N!");
    return this->M[v1-1][v2-1];
}

bool GraphContainer::GetElementAdjacencyMatrix(unsigned index) const
{
    if (index >= this->NTimesNMinusOneDiv2)
        throw std::invalid_argument("GetElementAdjacencyMatrix 0 <= index <= N(N-1)/2!");
    return this->GetElementAdjacencyMatrix(this->RowM[index], this->ColM[index]);
}

void GraphContainer::SetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("SetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    this->M[v1-1][v2-1] = true;
}

void GraphContainer::UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("UnsetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    this->M[v1-1][v2-1] = false;
}

int GraphContainer::GetVertexOrder(unsigned int v) const
{
    if (v > this->N || v < 1)
        throw std::invalid_argument("GetVertexOrder requires N >= v >= 1!");
    return this->VertexOrder[v-1];
}

void GraphContainer::SetVertexOrder(unsigned int v, int order)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("SetVertexOrder requires N >= v >= 1!");
    if (order > this->N || order < 0)
        throw std::invalid_argument("SetVertexOrder requires N >= order >= 1!");
    this->VertexOrder[v-1] = order;
}
