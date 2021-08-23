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

   GraphContainer(int n, int m); /// constructor with no graph (need a further call to set graph)

   GraphContainer(int n, int m, graph *g);

   void SetGraphFromNauty(graph *g); /// adjacency matrix, bond counts, and vertex order counts from dense nauty graph (assumes N and MWords set!)

   bool GetElementAdjacencyMatrix(unsigned int v1, unsigned v2) const;

   bool GetElementAdjacencyMatrix(unsigned index) const;

   int GetVertexOrder(unsigned int v) const; /// get the order of a given vertex v

   int GetN() const { return this->N; }

   int GetL() const { return this->L; }

   int GetNTimesNMinusOneDiv2() const { return this->NTimesNMinusOneDiv2; }

   int GetRowM(int index) const { return this->RowM[index]; } /// TODO: add checks

   int GetColM(int index) const { return this->ColM[index]; } /// TODO: add checks

};

#endif // GRAPHCONTAINER_H
