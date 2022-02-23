#ifndef WRAPPERROUTINESFOREMBEDDING_H
#define WRAPPERROUTINESFOREMBEDDING_H

#include "VertexEmbedList.h"
#include "GraphContainer.h"
#include "GraphEmbedder.h"
#include <iostream>

namespace WrapperRouintesForEmbedding {
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int> > ComputeRootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container, MaxInteractionLength correlatorLength);
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int> > ComputeUnrootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container);
std::vector<std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>>> ComputeRootedCanonicalEmbeddingsAndCountsCubic(const GraphContainer& container, MaxInteractionLength embeddingLength, MaxInteractionLength correlatorLength);
std::vector<std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>>> ComputeUnrootedCanonicalEmbeddingsAndCountsCubic(const GraphContainer& container, MaxInteractionLength embeddingLength);
}

#endif // WRAPPERROUTINESFOREMBEDDING_H
