#ifndef GRAPHGENERATOR_H
#define GRAPHGENERATOR_H

#include <iostream>

#include "GenericConnectedGraph.h"

class GraphParameters /// graph parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// max order of vertices

    unsigned int LMax; /// max number of bonds

    bool Disconnected; /// do we allow disconnected diagrams

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    GraphParameters(int argc, char *argv[]);

    unsigned int GetN() const { return this->N; }

    unsigned int GetLMax()  const { return this->LMax; }

    bool AllowDisconnected() const { return this->Disconnected; }

};

class GenericConnectedUndirectedGraphGenerator /// generate all generic connected undirected graphs
{

private:
    GraphParameters Parameters;

    GenericUndirectedConnectedGraph Graph;

    void AddBond(unsigned int v1, unsigned int v2) { this->Graph.AddBond(v1+1, v2+1); } /// v1, v2 start from ZERO!

    void RemoveBond(unsigned int v1, unsigned int v2) { this->Graph.RemoveBond(v1+1, v2+1); } /// v1, v2 start from ZERO!

    bool IsCanonical(int col, bool verbose=false) { return this->Graph.IsCanonical(col, verbose); } /// v2 starts from ZERO!

    bool IsConnected() { return this->Graph.IsConnected(); }

    unsigned int GetLMax() const { return this->Graph.GetLMax(); }

    unsigned int GetN() const { return this->Graph.GetN(); }

    int ComputeCurrentKey() { return this->Graph.ComputeCurrentKey(); }

    int GetRowForGenerate(int index) { return this->Graph.GetRowM(index)-1; } /// start from ZERO!

    int GetColForGenerate(int index) { return this->Graph.GetColM(index)-1; }  /// start from ZERO!

public:
    GenericConnectedUndirectedGraphGenerator(int argc, char *argv[]);

    void Generate(bool verbose=false); /// routine which does generation

    void GenerateNew(bool verbose=false); /// routine which does generation

};

#endif // GRAPHGENERATOR_H
