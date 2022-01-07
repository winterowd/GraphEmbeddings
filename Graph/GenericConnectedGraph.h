#ifndef GENERICCONNECTEDGRAPH_H
#define GENERICCONNECTEDGRAPH_H

#include <vector>

class GenericUndirectedConnectedGraph
{
private:
    /**** private variables *****/
    unsigned int N; // fixed order

    unsigned int LMax; // max number of bonds

    unsigned int NTimesNMinusOneDiv2; /// precompute N(N-1)2

    unsigned int CurrentNbrL; ///  current number of bonds

    unsigned int CanonicalKey; /// key for canonical labeling of a graph with L bonds and N vertices

    std::vector<int> VertexOrder; /// vector of vertex orders

    std::vector<int> VertexOrderTestCanonical; /// vector of vertex orders for canonical test

    std::vector<int> VertexOrderCounts; /// vector of counts of vertices with a given order 1 <= order <= N-1

    std::vector<int> VertexOrderCountsTestCanonical; /// vector of counts of vertices with a given order 1 <= order <= N-1 for canonical test

    std::vector<std::vector<bool>> M; /// NxN adjacency matrix

    std::vector<std::vector<bool>> MTestCanonical; /// copy of adjacency matrix for internal purposes

    /// debugging
    void PrintM();
    void PrintMTestCanonical();
    void PrintVertexOrders();
    void PrintVertexOrdersTestCanonical();

    /// precompute these in constructor
    std::vector<int> RowM;
    std::vector<int> ColM;

     /**** private methods *****/

    bool VertexOrdersConsistentWithAdjacencyMatrix(); /// check vertex orders agree with sums over rows

    bool AdjacencyMatrixOK(); /// zeros on diagonal and symmetric?

    bool VertexOrderCountsConsistent(); /// counts consistent?

    int ComputeVertexOrder(unsigned int v); /// compute the vertex order by summing over roq

    int ComputeVertexOrderTestCanonical(unsigned int v); // compute the vertex order of test canonical by summing over row

    int ComputeNbrBonds(); /// compute total number of bonds from adjacency matrix

    bool VertexOrderHistogramsEqual(); /// compare vertex order counts (histograms) for graph and test canonical

    bool NeighborVertexOrderHistogramsEqual(); /// compute the

    bool IsCanonicalPegsInHoles(int col, int key, int verbose=false);

    int ComputeKeyFromVector(const std::vector<bool>& vec); /// compute integer key using binary string provided by vec

    int ComputeKeyCol(int v2); /// compute integer key from

    unsigned int FindLastVertex() const;

    /***** private accessors ******/

    void SetVertexOrder(unsigned int v, int order);

    void SetVertexOrderTestCanonical(unsigned int v, int order);

    void IncrementVertexOrder(unsigned int v);

    void DecrementVertexOrder(unsigned int v);

    void SetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

    void SetElementAdjacenyMatrixTestCanonical(unsigned int v1, unsigned v2);

    void UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

    void UnsetElementAdjacenyMatrixTestCanonical(unsigned int v1, unsigned v2);

    /// decrease count for vertices of the same order as v
    void DecrementVertexOrderCount(unsigned int v);

    /// increase count for vertices of the same order as v
    void IncrementVertexOrderCount(unsigned int v);

    int ComputeTestCanonicalKey(); /// compute integer key using binary string from upper triangular part of matrix

    bool GetElementAdjacencyMatrixTestCanonical(unsigned int v1, unsigned v2) const;

public:

    GenericUndirectedConnectedGraph(unsigned int n, unsigned int lMax); /// constructor

    bool GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const;

    void SetGraphFromAdjacencyMatrix(const std::vector<std::vector<bool>>& m); /// set up graph from an adjacency matrix

    void AddBond(unsigned int v1, unsigned v2); /// add bond between v1 and v2

    void RemoveBond(unsigned int v1, unsigned v2); /// remove bond between v1 and v2

    int GetVertexOrder(unsigned int v) const;

    int GetVertexOrderTestCanonical(unsigned int v) const;

    int GetCurrentNbrL() const { return this->CurrentNbrL; }

    unsigned int GetLMax() const { return this->LMax; }

    unsigned int GetN() const { return this->N; }

    unsigned int GetNTimesNMinusOneDiv2() const { return this->NTimesNMinusOneDiv2; }

    int GetRowM(unsigned int index) const;  /// add checks for range of index

    int GetColM(unsigned int index) const;  /// add checks for range of index

    bool IsCanonical(int col, bool verbose=false); /// is the graph canonical?

    bool IsConnected(bool testCanonical=false); /// is the graph connected?

    bool IsEmbeddableCubic(); /// need to write this

    void ClearGraph(); /// clear the graph and all auxiliary variables

    int ComputeCurrentKey(); /// compute integer key using binary string from upper triangular part of matrix

    bool IsConsistent() { if (this->AdjacencyMatrixOK() && this->VertexOrdersConsistentWithAdjacencyMatrix()) return true; return false; } /// checks consistency of adjacency matrix and vertex orders

    friend std::ostream& operator<< (std::ostream& stream, const GenericUndirectedConnectedGraph& graph); /// output graph

    std::vector<std::vector<bool>> GetM() { return this->M; } /// return adjacency matrix by value for debugging

};

#endif // GENERICCONNECTEDGRAPH_H
