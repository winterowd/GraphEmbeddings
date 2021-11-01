#ifndef GRAPHCONTAINER_H
#define GRAPHCONTAINER_H

#include <vector>

extern "C" {
#include "nauty.h"
#include "gtools.h"
}

class GraphContainer
{

private:
    /**** private variables *****/
    int N; // fixed order

    int MWords; /// nauty result from SETWORDSNEEDED(n)

    int L; // number of bonds

    bool StoreRooted; /// flag to store rooted vertices

    int NbrRooted; /// number of rooted vertices

    std::vector<int> RootedVertices; /// container for storing labels of rooted vertices

    int NTimesNMinusOneDiv2; /// precompute N(N-1)2

    std::vector<int> VertexOrder; /// vector of vertex orders

    std::vector<std::vector<bool>> M; /// NxN adjacency matrix

    /// precompute these in constructor
    std::vector<int> RowM;
    std::vector<int> ColM;

    /**** private methods *****/

   bool VertexOrdersConsistentWithAdjacencyMatrix(); /// DEBUGGING: check vertex orders agree with sums over rows

   bool AdjacencyMatrixOK(); /// DEBUGGING: zeros on diagonal and symmetric?

   void SetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

   void UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

   int ComputeVertexOrder(unsigned int v);

   void SetRowAndColMapping();

   void SetVertexOrder(unsigned int v, int order);

public: /// later make this private?
   void PrintM() const;
   void PrintVertexOrders() const;

public:

   GraphContainer(int n, int m, bool storeRooted=false, int nbrRooted=1); /// constructor with no graph (need a further call to set graph)

   GraphContainer(int n, int m, graph *g, bool storeRooted=false, int nbrRooted=1);

   void SetGraphFromDenseNauty(graph *g); /// adjacency matrix, bond counts, and vertex order counts from dense nauty graph (assumes N and MWords set!)

   void GetDenseNautyFromGraph(graph *g) const; /// adjacency matrix to dense nauty graph

   void ColoredCanonicalRelabeling(int *labCanon, int v1, int v2=-1, bool verbose=false);

   void RelabelVertices(const std::vector<int>& newLabels);

   bool GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const;

   bool GetElementAdjacencyMatrix(unsigned index) const;

   int ComputeCurrentKey() const;

   bool StoringRooted() const { return this->StoreRooted; }

   int GetNbrRooted() const;

   int GetVertexOrder(unsigned int v) const; /// get the order of a given vertex v

   void SetRootedVertex(int index, int label);

   int GetRootedVertex(int index) const;

   int GetN() const { return this->N; }

   int GetL() const { return this->L; }

   int GetNTimesNMinusOneDiv2() const { return this->NTimesNMinusOneDiv2; }

   int GetRowM(int index) const;

   int GetColM(int index) const;

};

bool operator==(const GraphContainer& lhs, const GraphContainer& rhs); /// comparison operator
bool operator!=(const GraphContainer& lhs, const GraphContainer& rhs);
bool operator<(const GraphContainer& lhs, const GraphContainer& rhs);

#endif // GRAPHCONTAINER_H
