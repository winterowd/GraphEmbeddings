#ifndef WRAPPERROUTINESFOREMBEDDING_H
#define WRAPPERROUTINESFOREMBEDDING_H

#include "VertexEmbedList.h"
#include "GraphContainer.h"
#include "GraphEmbedder.h"
#include <iostream>

namespace WrapperRouintesForEmbedding {
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> ComputeRootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container, MaxInteractionLength correlatorLength);
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> ComputeUnrootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container);
std::pair<GraphContainer, std::vector<std::tuple<std::vector<VertexEmbedList>, std::vector<int>, std::vector<int>>>> ComputeRootedCanonicalEmbeddingsAndCountsCubic(const GraphContainer& container, MaxInteractionLength maxEmbeddingLength, MaxInteractionLength correlatorLength, int maxManhattanDistance);
std::pair<GraphContainer, std::vector<std::tuple<std::vector<VertexEmbedList>, std::vector<int>, std::vector<int>>>> ComputeUnrootedCanonicalEmbeddingsAndCountsCubic(const GraphContainer& container, MaxInteractionLength maxEmbeddingLength);
}

#endif // WRAPPERROUTINESFOREMBEDDING_H
