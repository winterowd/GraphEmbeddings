#include "WrapperRoutinesForEmbedding.h"

/// wrapper which creates GraphEmbedder object and returns the different embeddings for a two-rooted graph (\Lambda_1 or NN only!)
/// @param container: container for two-rooted graph
/// @param correlatorLength: separation of two-rooted vertices
/// @return tuple where the first argument is the symmetry factor, second is a vector of canoncial embeddings, vector of canonical weak embedding number
std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> WrapperRouintesForEmbedding::ComputeRootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container, MaxInteractionLength correlatorLength)
{
    if (container.GetNbrRooted()!=2)
        throw std::invalid_argument("ComputeRootedCanonicalEmbeddingsAndCountsCubicNN expects a two-rooted container!\n");

    for (int i=0; i<2; ++i)
        if (container.GetRootedVertex(i)!=i)
            throw std::invalid_argument("ComputeRootedCanonicalEmbeddingsAndCountsCubicNN rooted vertex i to have label i!\n");

    // prepare arguments for GraphEmbedder object
    std::vector<std::string> arguments;

    arguments.push_back("Foo"); /// dummy first argument
    arguments.push_back("-n"); /// number of vertices
    arguments.push_back(std::to_string(container.GetN()));
    arguments.push_back("-g"); /// g6 string
    arguments.push_back(container.GetG6String());
    arguments.push_back("-c"); /// correlator flag
    arguments.push_back("-d"); /// max embedding length
    arguments.push_back("NN");
    arguments.push_back("-m"); /// correlator length (NN)
    arguments.push_back(StringFromMaxInteractionLength(correlatorLength));
    arguments.push_back("-l"); /// cubic lattice
    arguments.push_back("Cubic");

    std::vector<char*> argv;
    for (const auto& arg: arguments)
        argv.push_back((char*)arg.data());
    argv.push_back(nullptr);

    GraphEmbedder embedder(argv.size()-1,argv.data());

    return embedder.ComputeCanonicalEmbeddingsAndCountsNN(container);
}

std::tuple<GraphContainer, std::vector<VertexEmbedList>, std::vector<int>> WrapperRouintesForEmbedding::ComputeUnrootedCanonicalEmbeddingsAndCountsCubicNN(const GraphContainer& container)
{
    if (container.GetNbrRooted()!=0)
        throw std::invalid_argument("ComputeUnrootedCanonicalEmbeddingsAndCountsCubicNN expects an unrooted container!\n");

    // prepare arguments for GraphEmbedder object
    std::vector<std::string> arguments;

    arguments.push_back("Foo"); /// dummy first argument
    arguments.push_back("-n"); /// number of vertices
    arguments.push_back(std::to_string(container.GetN()));
    arguments.push_back("-g"); /// g6 string
    arguments.push_back(container.GetG6String());
    arguments.push_back("-d"); /// max embedding length
    arguments.push_back("NN");
    arguments.push_back("-l"); /// cubic lattice
    arguments.push_back("Cubic");

    std::vector<char*> argv;
    for (const auto& arg: arguments)
        argv.push_back((char*)arg.data());
    argv.push_back(nullptr);

    GraphEmbedder embedder(argv.size()-1,argv.data());

    return embedder.ComputeCanonicalEmbeddingsAndCountsNN(container);
}
