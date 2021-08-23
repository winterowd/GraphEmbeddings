#include "GenericConnectedGraph.h"

#include <tuple>
#include <iostream>

void GenericUndirectedConnectedGraph::SetGraphFromAdjacencyMatrix(const std::vector<std::vector<bool>>& m)
{
    if (m.size() != this->M.size())
        throw std::invalid_argument("SetGraphFromAdjacencyMatrix M not of the proper size!");

    for (int i=0; i<m.size(); ++i)
        if (m[i].size()!=this->M[i].size())
            throw std::invalid_argument("SetGraphFromAdjacencyMatrix M not of the proper size!");

    this->M = m; /// copy adjacency matrix

    this->CurrentNbrL = this->ComputeNbrBonds(); /// update the number of bonds

    /// perform checks on number of bonds
    if (this->CurrentNbrL>this->LMax)
        throw std::invalid_argument("SetGraphFromAdjacencyMatrix number of bonds exceeds LMax!");

    /// set the vertex orders
    for (unsigned int v=1; v<=this->N; ++v)
        this->SetVertexOrder(v, this->ComputeVertexOrder(v));

    /// set the vertex counts
    for (unsigned int v=1; v<=this->N; ++v)
        this->VertexOrderCounts[this->VertexOrder[v-1]]++;

    if (!this->IsConnected())
        throw std::invalid_argument("SetGraphFromAdjacencyMatrix graph is disconnected!");

}

int GenericUndirectedConnectedGraph::GetRowM(unsigned int index) const
{
    if (index < 0 || index > this->RowM.size())
        throw std::invalid_argument("GetRowM index out of range!");
    return this->RowM[index];
}

int GenericUndirectedConnectedGraph::GetColM(unsigned int index) const
{
    if (index < 0 || index > this->RowM.size())
        throw std::invalid_argument("GetColM index out of range!");
    return this->ColM[index];
}

void GenericUndirectedConnectedGraph::ClearGraph()
{
    for (int i=0; i<M.size(); ++i)
        for (int j=0; j<M.size(); ++j)
            M[i][j] = false;

    for (int i=0; i<VertexOrder.size(); ++i)
        VertexOrder[i] = 0;

    VertexOrderCounts[0] = this->N;
    for (int i=1; i<VertexOrderCounts.size(); ++i)
        VertexOrderCounts[i] = 0;

    this->CurrentNbrL = 0;
    this->CanonicalKey = 0;
}

bool GenericUndirectedConnectedGraph::GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("GetElementAdjacencyMatrix 1 <= v1, v2 <= N!");
    return this->M[v1-1][v2-1];
}

bool GenericUndirectedConnectedGraph::GetElementAdjacencyMatrixTestCanonical(unsigned int v1, unsigned v2) const
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("GetElementAdjacencyMatrixTestCanonical 1 <= v1, v2 <= N!");
    return this->MTestCanonical[v1-1][v2-1];
}

void GenericUndirectedConnectedGraph::SetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("SetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    /// warning if it is already occupied
    if (this->GetElementAdjacencyMatrix(v1,v2))
        std::cerr << "WARNING: In SetElementAdjacenyMatrix setting a link that was already present!" << std::endl;
    this->M[v1-1][v2-1] = true;
}

void GenericUndirectedConnectedGraph::SetElementAdjacenyMatrixTestCanonical(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("SetElementAdjacenyMatrixTestCanonical 1 <= v1, v2 <= N!");
    this->MTestCanonical[v1-1][v2-1] = true;
}

void GenericUndirectedConnectedGraph::UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("UnsetElementAdjacenyMatrix 1 <= v1, v2 <= N!");
    /// warning if it is already vacant
    if (!this->GetElementAdjacencyMatrix(v1,v2))
        std::cerr << "WARNING: In UnsetElementAdjacenyMatrix unsetting a link that was unset!" << std::endl;
    this->M[v1-1][v2-1] = false;
}

void GenericUndirectedConnectedGraph::UnsetElementAdjacenyMatrixTestCanonical(unsigned int v1, unsigned v2)
{
    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("UnsetElementAdjacenyMatrixTestCanonical 1 <= v1, v2 <= N!");
    this->MTestCanonical[v1-1][v2-1] = false;
}

void GenericUndirectedConnectedGraph::DecrementVertexOrderCount(unsigned int v)
{
    this->VertexOrderCounts[this->VertexOrder[v-1]]--;
}

void GenericUndirectedConnectedGraph::IncrementVertexOrderCount(unsigned int v)
{
    this->VertexOrderCounts[this->VertexOrder[v-1]]++;
}

void GenericUndirectedConnectedGraph::AddBond(unsigned int v1, unsigned v2)
{
    if (v1 == v2)
        throw std::invalid_argument("AddBond must connect two different vertices v1 and v2!");

    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("AddBond requires N >= v1, v2 >= 1!");

    /// set elements of adjacency matrix
    this->SetElementAdjacenyMatrix(v1,v2);
    this->SetElementAdjacenyMatrix(v2,v1);

    /// update vertex order counts for v1 and v2's old order
    this->DecrementVertexOrderCount(v1);
    this->DecrementVertexOrderCount(v2);

    /// update vertex order for v1 and v2
    this->IncrementVertexOrder(v1);
    this->IncrementVertexOrder(v2);

    /// update vertex order counts for v1 and v2's new order
    this->IncrementVertexOrderCount(v1);
    this->IncrementVertexOrderCount(v2);

    /// update current number of links
    this->CurrentNbrL++;

    if (!(this->VertexOrdersConsistentWithAdjacencyMatrix() && this->AdjacencyMatrixOK() && this->VertexOrderCountsConsistent()))
        std::cerr << "Graph inconsistent after AddBond!\n";
}

void GenericUndirectedConnectedGraph::RemoveBond(unsigned int v1, unsigned v2)
{
    if (v1 == v2)
        throw std::invalid_argument("RemoveBond must connect two different vertices v1 and v2!");

    if (v1 > this->N || v1 < 1 || v2 > this->N || v2 < 1)
        throw std::invalid_argument("RemoveBond requires N >= v1, v2 >= 1!");

    /// set elements of adjacency matrix
    this->UnsetElementAdjacenyMatrix(v1,v2);
    this->UnsetElementAdjacenyMatrix(v2,v1);

    /// update vertex order counts for v1 and v2's old order
    this->DecrementVertexOrderCount(v1);
    this->DecrementVertexOrderCount(v2);

    /// update vertex order for v1 and v2
    this->DecrementVertexOrder(v1);
    this->DecrementVertexOrder(v2);

    /// update vertex order counts for v1 and v2's new order
    this->IncrementVertexOrderCount(v1);
    this->IncrementVertexOrderCount(v2);

    /// update current number of links
    this->CurrentNbrL--;

    if (!(this->VertexOrdersConsistentWithAdjacencyMatrix() && this->AdjacencyMatrixOK() && this->VertexOrderCountsConsistent()))
        std::cerr << "Graph inconsisten after RemoveBond!\n";
}

int GenericUndirectedConnectedGraph::GetVertexOrder(unsigned int v) const
{
    if (v > this->N || v < 1)
        throw std::invalid_argument("GetVertexOrder requires N >= v >= 1!");
    return this->VertexOrder[v-1];
}

int GenericUndirectedConnectedGraph::GetVertexOrderTestCanonical(unsigned int v) const
{
    if (v > this->N || v < 1)
        throw std::invalid_argument("GetVertexOrderTestCanonical requires N >= v >= 1!");
    return this->VertexOrderTestCanonical[v-1];
}

void GenericUndirectedConnectedGraph::SetVertexOrder(unsigned int v, int order)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("SetVertexOrder requires N >= v >= 1!");
    if (order > this->N || order < 0)
        throw std::invalid_argument("SetVertexOrder requires N >= order >= 1!");
    this->VertexOrder[v-1] = order;
}

void GenericUndirectedConnectedGraph::SetVertexOrderTestCanonical(unsigned int v, int order)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("SetVertexOrderTestCanonical requires N >= v >= 1!");
    if (order > this->N || order < 0)
        throw std::invalid_argument("SetVertexOrderTestCanonical requires N >= order >= 1!");
    this->VertexOrderTestCanonical[v-1] = order;
}

void GenericUndirectedConnectedGraph::IncrementVertexOrder(unsigned int v)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("IncrementVertexOrder requires N >= v >= 1!");
    if (this->VertexOrder[v-1] == this->N)
        throw std::invalid_argument("IncrementVertexOrder cannot increment vertex order to be greater than N!");
    this->VertexOrder[v-1]++;
}

void GenericUndirectedConnectedGraph::DecrementVertexOrder(unsigned int v)
{
    if (v>this->N || v<=0)
        throw std::invalid_argument("DecrementVertexOrder requires N >= v >= 1!");
    if (this->VertexOrder[v-1] == 0)
        throw std::invalid_argument("DecrementVertexOrder cannot increment vertex order to be negative!");
    this->VertexOrder[v-1]--;
}

unsigned int GenericUndirectedConnectedGraph::FindLastVertex() const
{
    unsigned int vLast = 0;
    for (int i=this->NTimesNMinusOneDiv2-1; i>=0; --i)
    {
        if (this->GetElementAdjacencyMatrix(this->RowM[i],this->ColM[i]))
        {
            vLast = this->ColM[i];
            break;
        }
    }
    return vLast;
}

/// depth first search for connectivity
bool GenericUndirectedConnectedGraph::IsConnected(bool testCanonical)
{
    unsigned int vLast = this->FindLastVertex();

    std::vector<int> q;
    std::vector<bool> visited(vLast, false);

    visited[0] = true;
    q.push_back(0);
    while (!q.empty())
    {
        int j = q.back();
        q.pop_back();
        for (int k=0; k<this->N; ++k)
        {
            if (j!=k)
            {
                if (testCanonical)
                {
                    if (this->GetElementAdjacencyMatrixTestCanonical(j+1,k+1) && !visited[k])
                    {
                        visited[k] = true;
                        q.push_back(k);
                    }
                }
                else
                {
                    if (this->GetElementAdjacencyMatrixTestCanonical(j+1,k+1) && !visited[k])
                    {
                        visited[k] = true;
                        q.push_back(k);
                    }
                }
            }
        }
    }
    for (auto v: visited)
        if (!v)
            return false;
   return true;
}


bool GenericUndirectedConnectedGraph::VertexOrderHistogramsEqual()
{
    for (int i=0; i<this->VertexOrderCounts.size(); ++i)
    {
        if (this->VertexOrderCounts[i] != this->VertexOrderCountsTestCanonical[i])
            return false;
    }
    return true;
}

bool GenericUndirectedConnectedGraph::NeighborVertexOrderHistogramsEqual()
{
    /// for each possible vertex order, store counts for vertex orders of neighbors
    std::vector<std::vector<int>> histogram(this->N-1, std::vector<int>(this->N-1,0));
    std::vector<std::vector<int>> histogramTest(this->N-1, std::vector<int>(this->N-1,0));

    for (int v1=1; v1<=this->N; ++v1)
    {
        for (int v2=v1+1; v2<=this->N; ++v2)
        {
            if (this->GetElementAdjacencyMatrix(v1,v2))
            {
                histogram[this->GetVertexOrder(v1)-1][this->GetVertexOrder(v2)-1]++; ///
                histogram[this->GetVertexOrder(v2)-1][this->GetVertexOrder(v1)-1]++;
            }

            if (this->GetElementAdjacencyMatrixTestCanonical(v1,v2))
            {
                histogramTest[this->GetVertexOrderTestCanonical(v1)-1][this->GetVertexOrderTestCanonical(v2)-1]++; ///
                histogramTest[this->GetVertexOrderTestCanonical(v2)-1][this->GetVertexOrderTestCanonical(v1)-1]++;
            }
        }
    }

    /*for (int i=0; i<histogram.size(); ++i)
        for (int j=0; j<histogram[i].size(); ++j)
            std::cout << "VertexOrder " << i+1 << " has " << histogram[i][j] << " counts of neighbors of order " << j+1 <<  " vs " << histogramTest[i][j] << "\n";*/

    for (int i=0; i<histogram.size(); ++i)
        for (int j=0; j<histogram[i].size(); ++j)
            if (histogram[i][j] != histogramTest[i][j])
                return false;

    return true;
}

bool GenericUndirectedConnectedGraph::IsCanonicalPegsInHoles(int col, int key, int verbose)
{
    unsigned int m = this->GetCurrentNbrL();
    unsigned int n = this->NTimesNMinusOneDiv2;

    if (m>n)
        throw std::invalid_argument("IsCanonicalPegsInHoles: M must be less than or equal to N!");

    /// pegs in holes with 1s and 0s
    std::vector<int> s(m, -1);
    std::vector<bool> occ(n, false);

    int count = 0;
    int sk;
    int k=0;
    while (1)
    {

        if (s[k] < 0 && k > 0)
            sk = s[k-1]+1;
        else
            sk = s[k]+1; /// next location of kth peg
        s[k] = sk; /// set kth peg

        if (sk>=n) /// check if we went over the number of holes or it is disconnected
        {
            s[k] = -1; /// reset kth peg
            k--; /// go back to previous peg (k-1)
            if (k==0) /// we are at the end if the first peg has reached past the last hole
                break;
            occ[s[k]] = false; /// mark as unoccupied as we intend to move (k-1) peg to new location
        }
        else
        {
            if (!occ[sk]) /// check if occupied
            {
                occ[sk] = true; /// set occupied if it is
                k++; /// go to next peg
            }
            else
                continue; /// if it is occupied, move kth peg again

            if (k>=m) /// have we successfully placed all the pegs?
            {
                count++;
                if (verbose)
                {
                    std::cout << "FOUND A CONFIG! " << count << std::endl;
                    for (int i=0; i<m; ++i)
                        std::cout << s[i] << " ";
                    std::cout << "\n";
                    for (int i=0; i<n; ++i)
                        std::cout << occ[i] << " ";
                    std::cout << "\n";
                }

                /// update MTestCanonical
                for (int i=0; i<occ.size(); ++i)
                {
                    if (occ[i])
                    {
                        this->SetElementAdjacenyMatrixTestCanonical(RowM[i], ColM[i]);
                        this->SetElementAdjacenyMatrixTestCanonical(ColM[i], RowM[i]);
                    }
                    else
                    {
                        this->UnsetElementAdjacenyMatrixTestCanonical(RowM[i], ColM[i]);
                        this->UnsetElementAdjacenyMatrixTestCanonical(ColM[i], RowM[i]);
                    }
                }

                /// clear histogram for test canonical
                for (auto i=0; i<this->VertexOrderCountsTestCanonical.size(); ++i)
                    this->VertexOrderCountsTestCanonical[i] = 0;

                /// set the vertex orders and counts
                for (int v=1; v<=this->N; ++v)
                {
                    this->SetVertexOrderTestCanonical(v, this->ComputeVertexOrderTestCanonical(v));
                    this->VertexOrderCountsTestCanonical[this->VertexOrderTestCanonical[v-1]]++;
                }

                /// compute key associated to occ and compare with key from original column
                /// check that vertex order histograms are equal
                auto testCanonicalKey = this->ComputeTestCanonicalKey();
                if (this->ComputeTestCanonicalKey() > key)
                {
                    if (verbose)
                        std::cout << "FOUND GREATER KEY! " << testCanonicalKey << " vs " << key << "\n";
                    if (this->VertexOrderHistogramsEqual() && this->IsConnected(true))
                    {
                        if (verbose)
                            std::cout << "VERTEX_ORDER_COUNTS_EQUAL! " << this->CurrentNbrL << "\n";

                        if (this->NeighborVertexOrderHistogramsEqual())
                        {
                            if (verbose)
                                std::cout << "NEIGHBOR VERTEX_ORDER_COUNTS_EQUAL!\n";
                            return false;
                        }
                        else
                            if (verbose)
                                std::cout << "NEIGHBOR VERTEX_ORDER_COUNTS_NOT_EQUAL!\n";
                    }
                    else
                        if (verbose)
                            std::cout << "VERTEX_ORDER_COUNTS_NOT_EQUAL OR DISCONNECTED!\n";
                }
                else
                    if (verbose)
                        std::cout << "NEW KEY IS LESS THAN OR EQUAL TO CURRENT KEY! " << testCanonicalKey << " vs " << key << "\n";

                k--; /// go back (we are over the size of the array s)
                occ[s[k]] = false; /// mark as unoccupied as we intend to move kth peg to new location
            }

        }
    }
    return true;
}


/// a new bond added to column v2 of adjacency matrix: 1 <= v2 <= N
/// assumes previous graph was canonical
bool GenericUndirectedConnectedGraph::IsCanonical(int col, bool verbose)
{

    /// copy adjacency matrix
    this->MTestCanonical = this->M;

    /// pegs in holes with 1s and 0s
    if (this->IsCanonicalPegsInHoles(col, this->ComputeCurrentKey(), verbose))
    {
        if (verbose)
            std::cout << "Graph is canonical!" << std::endl;
        return true;
    }
    else
    {
        if (verbose)
            std::cout << "Graph is not canonical!" << std::endl;
        return false;
    }
}

int GenericUndirectedConnectedGraph::ComputeVertexOrder(unsigned int v)
{
    int result = 0;
    for (int i=0; i<this->M[v-1].size(); ++i) /// sum along row of adjacency matrix
        result += this->M[v-1][i];
    return result;
}

int GenericUndirectedConnectedGraph::ComputeVertexOrderTestCanonical(unsigned int v)
{
    int result = 0;
    for (int i=0; i<this->M[v-1].size(); ++i) /// sum along row of adjacency matrix
        result += this->MTestCanonical[v-1][i];
    return result;
}

int GenericUndirectedConnectedGraph::ComputeNbrBonds()
{
    int nbrBonds = 0;
    for (int i=0; i<this->M.size(); ++i)
        for (int j=i+1; j<this->M.size(); ++j)
            if (this->M[i][j])
                nbrBonds++;
    return nbrBonds;
}

int GenericUndirectedConnectedGraph::ComputeKeyFromVector(const std::vector<bool>& vec)
{
    int result = 0;
    for (int i=0; i<vec.size(); ++i)
        result += vec[i] << i;
    return result;
}

bool GenericUndirectedConnectedGraph::VertexOrdersConsistentWithAdjacencyMatrix()
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

bool GenericUndirectedConnectedGraph::AdjacencyMatrixOK()
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

bool GenericUndirectedConnectedGraph::VertexOrderCountsConsistent()
{
    auto result = true;
    std::vector<int> tempCounts(this->N, 0);
    for (int v=1; v<=this->N; ++v)
        tempCounts[this->GetVertexOrder(v)]++;
    for (int i=0; i<tempCounts.size(); ++i)
        if (tempCounts[i]!=this->VertexOrderCounts[i])
        {
            std::cout << "ERROR: VertexOrderCountsConsistent counts  at " << i << " non consistent! " << tempCounts[i] << " vs " << this->VertexOrderCounts[i] << "\n";
            result = false;
        }
    return result;
}

int GenericUndirectedConnectedGraph::ComputeKeyCol(int v2)
{
    int result = 0;
    int count = 0;
    for (int v1=v2-1; v1>0; --v1)
    {
        result += this->GetElementAdjacencyMatrix(v1,v2) << count;
        count++;
    }
    return result;
}

int GenericUndirectedConnectedGraph::ComputeCurrentKey()
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

int GenericUndirectedConnectedGraph::ComputeTestCanonicalKey()
{
    int result = 0;
    int count = 0;
    for (int i=this->NTimesNMinusOneDiv2-1; i>=0; --i)
    {
       result += this->GetElementAdjacencyMatrixTestCanonical(this->RowM[i],this->ColM[i]) << count;
       count++;
    }
    return result;
}

void GenericUndirectedConnectedGraph::PrintM()
{
    for (int i=0; i<this->NTimesNMinusOneDiv2; ++i)
        std::cout << "M: " << i << " " << RowM[i] << " " << ColM[i] << " " << this->GetElementAdjacencyMatrix(RowM[i], ColM[i]) << "\n";
}

void GenericUndirectedConnectedGraph::PrintMTestCanonical()
{
    for (int i=0; i<this->NTimesNMinusOneDiv2; ++i)
        std::cout << "MTestCanonical: " << i << " " << RowM[i] << " " << ColM[i] << " " << this->GetElementAdjacencyMatrixTestCanonical(RowM[i], ColM[i]) << "\n";
}

void GenericUndirectedConnectedGraph::PrintVertexOrders()
{
    for (int i=1; i<=this->N; ++i)
        std::cout << "VERTEX_ORDER: " << i << " " << this->GetVertexOrder(i) << "\n";
}

void GenericUndirectedConnectedGraph::PrintVertexOrdersTestCanonical()
{
    for (int i=1; i<=this->N; ++i)
        std::cout << "VERTEX_ORDER_TEST_CANONICAL: " << i << " " << this->GetVertexOrderTestCanonical(i) << "\n";
}


GenericUndirectedConnectedGraph::GenericUndirectedConnectedGraph(unsigned int n, unsigned int lMax) :
    N(n),
    LMax(lMax),
    NTimesNMinusOneDiv2(n*(n-1)/2),
    CurrentNbrL(0),
    CanonicalKey(0),
    VertexOrder(n,0),
    VertexOrderTestCanonical(n,0),
    VertexOrderCounts(n,0),
    VertexOrderCountsTestCanonical(n,0),
    M(n, std::vector<bool>(n,false)),
    MTestCanonical(n, std::vector<bool>(n,false)),
    RowM(n*(n-1)/2),
    ColM(n*(n-1)/2)
{
    if (this->LMax<this->N-1 || this->LMax>this->NTimesNMinusOneDiv2)
        throw std::invalid_argument("GenericUndirectedConnectedGraph requires N(N-1)/2 >= lMax >= n-1!");

    VertexOrderCounts[0] = this->N;

    /// fill look up tables for upper triangular elements of adjacency matrix
    int b=-1;
    for (int i=1; i<n; ++i)
    {
        for (int j=0; j<i; ++j)
        {
            b++;
            this->RowM[b] = j+1;
            this->ColM[b] = i+1;
        }
    }
}

std::ostream& operator<< (std::ostream& stream, const GenericUndirectedConnectedGraph& graph)
{
    auto n = graph.FindLastVertex();
    stream << n << " " << graph.GetCurrentNbrL() << " ";
    for (int v1=1; v1<n; ++v1)
        for (int v2=v1+1; v2<=n; ++v2)
            if (graph.GetElementAdjacencyMatrix(v1,v2))
                stream << v1 << " " << v2 << " ";
    return stream;
}
