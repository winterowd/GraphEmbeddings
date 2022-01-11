#ifndef AUXILIARYROUTINESFORNAUTY_H
#define AUXILIARYROUTINESFORNAUTY_H

#include <vector>
#include <iostream>

#include "GraphContainer.h"

extern "C" {
#include "nausparse.h"
}

namespace AuxiliaryRoutinesForNauty {
void SetColoredPartition(int* c, int* lab, int* ptn, int n);
void SetColoredPartition(int n, int k, int* c, int* lab, int* ptn);
void SetColoredPartition(sparsegraph* g, int k, int* c, int* lab, int* ptn);
void SetVertexColors(int *c, const std::vector<int>& rootedVertices, int n, bool verbose=false);
GraphContainer GetCanonicalColoredGraphNauty(int n, const std::string& g6String, const std::vector<int>& rootedVertices); /// returns canonical colored graph
GraphContainer GetCanonicalGraphNauty(int n, const std::string& g6String); /// returns canonical graph
GraphContainer GetCanonicalGraphNauty(int n, graph *g); /// returns canonical graph
}

#endif // AUXILIARYROUTINESFORNAUTY_H
