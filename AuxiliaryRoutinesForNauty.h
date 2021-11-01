#ifndef AUXILIARYROUTINESFORNAUTY_H
#define AUXILIARYROUTINESFORNAUTY_H

#include <vector>
#include <iostream>

extern "C" {
#include "nausparse.h"
}

void SetColoredPartition(int* c, int* lab, int* ptn, int n);
void SetColoredPartition(int n, int k, int* c, int* lab, int* ptn);
void SetColoredPartition(sparsegraph* g, int k, int* c, int* lab, int* ptn);
void SetVertexColors(int *c, const std::vector<int>& rootedVertices, int n, bool verbose=false);

#endif // AUXILIARYROUTINESFORNAUTY_H
