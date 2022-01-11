#include "CanonicalGraphManager.h"

CanonicalGraphManager::CanonicalGraphManager(int lMax, bool rooted) :
    LMax(lMax),
    Rooted(rooted),
    TotalGraphs(0),
    CanonicalGraphs(this->LMax, std::vector<GraphContainer>())
{
    this->GenerateCanonicalGraphs();
}

/// run the generation of connected graphs of a fixed number of vertices and bonds using GraphGeneratorNauty
/// @param n: number of vertices
/// @param l: number of bonds
void CanonicalGraphManager::RunGraphGeneratorNauty(int n, int l)
{
    // prepare arguments for GraphGeneratorNauty object
    std::vector<std::string> arguments;

    arguments.push_back("-n "+std::to_string(n)); /// order of graph
    arguments.push_back("-l "+std::to_string(l)); /// number of bonds
    if (this->Rooted) /// rooted?
        arguments.push_back("-r "+std::to_string(1));

    std::vector<char*> argv;
    for (const auto& arg: arguments)
        argv.push_back((char*)arg.data());
    argv.push_back(nullptr);

    GraphGeneratorNauty generator(argv.size()-1,argv.data());

    generator.Generate();

    this->ImportFromFileAndCanonicalize(n, l, generator.GetOutputFilenameFixedOrder(this->Rooted));

}

///
void CanonicalGraphManager::ImportFromFileAndCanonicalize(int n, int l, std::string filename)
{

}

void CanonicalGraphManager::AddCanonicalGraph(const GraphContainer& g)
{
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    this->CanonicalGraphs[g.GetL()-1].push_back(g);
    this->TotalGraphs++;
}
