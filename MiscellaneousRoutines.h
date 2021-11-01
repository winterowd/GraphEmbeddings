#ifndef MISCELLANEOUSROUTINES_H
#define MISCELLANEOUSROUTINES_H

#include <vector>
#include <iostream>

#include "SquareLattice.h"
#include "VertexEmbedList.h"

void PegsInHoles(const unsigned int n, const unsigned int m, bool identicalPegs, bool verbose=false);
void GenerateAllPermutationsWithRepeats(std::vector<std::vector<int>> &lists, const std::string& s, std::vector<int>& pos, int n, const int& size);
bool DoesNotDoubleBack(const std::vector<int> &path);
int OppositeDir(const int& dir);

struct Site {
    int x, y, z;
};

bool operator==(const Site& lhs, const Site& rhs);
bool operator!=(const Site& lhs, const Site& rhs);
std::ostream& operator<< (std::ostream& stream, const Site& site);
Site MoveDir(const Site& start, int dir, int L);
bool IsBridgeValid(const Site& start, const Site& end, const std::vector<int>& bridgeDirs, std::vector<Site>& Bridge, int L, bool verbose=false);

void TestLexicographicalOrderingVertexEmbedList();

#endif // MISCELLANEOUSROUTINES_H
