#ifndef GRAPHCONTAINER_H
#define GRAPHCONTAINER_H

#include <vector>
#include <iostream>

extern "C" {
#include "nauty.h"
#include "gtools.h"
}

struct SiteCount {
    SiteCount(int in=0, int out=0) : NbrIn(in), NbrOut(out) {}
    int NbrIn;
    int NbrOut;
};

inline std::ostream& operator<<(std::ostream& os, const SiteCount& c)
{
    os << "(" << c.NbrIn << ", " << c.NbrOut << ")";
    return os;
}

struct UndirectedEdge {
    UndirectedEdge(int v1=0, int v2=0) : FirstVertex(v1), SecondVertex(v2) {}
    int FirstVertex;
    int SecondVertex;
};

inline std::ostream& operator<<(std::ostream& os, const UndirectedEdge& e)
{
    os << "(" << e.FirstVertex << ", " << e.SecondVertex << ")";
    return os;
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
    std::vector<int> RowM; /// starts at 1
    std::vector<int> ColM; /// starts at 1

    std::vector<UndirectedEdge> Edges; /// edges of graph

    std::string G6String; /// save this when setting g from dense nauty

    /// NOTE: only use this currently in GraphEmbedder::GetCanonicalGraphsAndCounts
    int SymmFactor; /// symmetry factor which can be saved for embedding

    /**** private methods *****/

   bool VertexOrdersConsistentWithAdjacencyMatrix(); /// DEBUGGING: check vertex orders agree with sums over rows

   bool AdjacencyMatrixOK(); /// DEBUGGING: zeros on diagonal and symmetric?

   void SetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

   void UnsetElementAdjacenyMatrix(unsigned int v1, unsigned v2);

   int ComputeVertexOrder(unsigned int v);

   void SetRowAndColMapping();

   void SetVertexOrder(unsigned int v, int order);

   int GetSizeFromG6(const std::string& g6string);

   void SetG6StringFromDenseNauty(graph *g);

   void SetGraphFromDenseNauty(graph *g); /// adjacency matrix, bond counts, and vertex order counts from dense nauty graph (assumes N and MWords set!)

   void ResetEdges();

public: /// later make this private?

   void PrintVertexOrders() const;

public:

   GraphContainer(int n, int m, graph *g, bool storeRooted=false, int nbrRooted=1);

   GraphContainer(int n, int m, const std::string& g6String, bool storeRooted=false, int nbrRooted=1);

   GraphContainer(int n, int m, const std::vector<UndirectedEdge>& edges, bool storeRooted=false, int nbrRooted=1);

   void GetDenseNautyFromGraph(graph *g) const; /// adjacency matrix to dense nauty graph

   void ColoredCanonicalRelabeling(int *labCanon, int v1, int v2=-1);

   void CanonicalRelabeling(int *labCanon);

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

   int GetSymmFactor() const { return this->SymmFactor; };

   std::string GetG6String() const { return this->G6String; }

   void SetSymmFactor(int symmFactor); /// only use this currently in GraphEmbedder::GetCanonicalGraphsAndCounts

   bool IsConnected(int vertexStart=0);

   UndirectedEdge GetEdge(int index);

   std::vector<UndirectedEdge> GetAllEdges();

   friend std::ostream& operator<< (std::ostream& stream, const GraphContainer& can);

};

bool operator==(const GraphContainer& lhs, const GraphContainer& rhs); /// comparison operator
bool operator!=(const GraphContainer& lhs, const GraphContainer& rhs);
bool operator<(const GraphContainer& lhs, const GraphContainer& rhs);

#endif // GRAPHCONTAINER_H
