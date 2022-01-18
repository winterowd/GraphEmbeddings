#include "CanonicalGraphManager.h"

CanonicalGraphManager::CanonicalGraphManager(int lMax, bool rooted) :
    LMax(lMax),
    NbrGraphs(0),
    NbrRootedGraphs(0),
    CanonicalGraphs(this->LMax, std::vector<GraphContainer>()),
    CanonicalRootedGraphs(this->LMax, std::vector<GraphContainer>())
{
    this->GenerateCanonicalGraphs();
}

void CanonicalGraphManager::GenerateCanonicalGraphs()
{
    for (int n=2; n<=this->LMax+1; ++n) ///
    {
        int connUpper = n*(n-1)/2;
        int minBondsFixedN = n-1;
        int maxBondsFixedN = (connUpper < this->LMax) ? connUpper : this->LMax;
        std::cout << "DEBUG_GENERATE N: " << n << " L: [" << minBondsFixedN << "," << maxBondsFixedN << "]\n";
        for (int l=minBondsFixedN; l<=maxBondsFixedN; ++l)
            this->GetGraphsNauty(n, l);
    }
}

/// run the generation of connected graphs of a fixed number of vertices and bonds using GraphGeneratorNauty
/// @param n: number of vertices
/// @param l: number of bonds
void CanonicalGraphManager::GetGraphsNauty(int n, int l)
{
    // prepare arguments for GraphGeneratorNauty object (unrooted)
    std::vector<std::string> arguments;

    arguments.push_back("Foo");
    arguments.push_back("-n"); /// order of graph
    arguments.push_back(std::to_string(n));
    arguments.push_back("-l"); /// number of bonds
    arguments.push_back(std::to_string(l));

    std::vector<char*> argv;
    for (const auto& arg: arguments)
        argv.push_back((char*)arg.data());
    argv.push_back(nullptr);

    GraphGeneratorNauty generator(argv.size()-1,argv.data());

    generator.Generate(); /// generate unrooted (n,l)

    /// add arguments for rooted graphs
    arguments.push_back("-r");
    arguments.push_back(std::to_string(1));

    std::vector<char*> argvRooted;
    for (const auto& arg: arguments)
        argvRooted.push_back((char*)arg.data());
    argvRooted.push_back(nullptr);

    GraphGeneratorNauty generatorRooted(argvRooted.size()-1,argvRooted.data());

    generatorRooted.Generate(); /// generate rooted (n,l)

    this->ImportFromFileAndCanonicalize(n, l, generator.GetOutputFilenameFixedOrder(), generator.GetOutputFilenameFixedOrder(true));

}

/// read in g6 from file: all graphs (rooted and unrooted) have n vertices and l bonds
void CanonicalGraphManager::ImportFromFileAndCanonicalize(int n, int l, std::string filenameUnrooted, std::string filenameRooted)
{
    std::ifstream unrootedFile(filenameUnrooted, std::ifstream::in); /// unrooted file
    if (!unrootedFile.is_open())
        throw std::logic_error("ImportFromFileAndCanonicalize: unrootedFile not open!\n");
    std::string tempLine;
    int m = SETWORDSNEEDED(n);
    int count=0;
    while (std::getline(unrootedFile, tempLine))
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)
#ifdef DEBUG
        auto temp = AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(n, g6String);
        if (temp.GetL()!=l)
            throw std::invalid_argument("ImportFromFileAndCanonicalize requires the correct number of bonds!\n");
        this->AddCanonicalGraph(temp);
#else
        this->AddCanonicalGraph(AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(n, g6String));
#endif
        count++;
    }
    unrootedFile.close();
#ifdef DEBUG
    std::cout << "IMPORTED_UNROOTED: (" << n << "," << l << "): " << count << "\n";
#endif


    std::ifstream rootedFile(filenameRooted, std::ifstream::in); /// unrooted file
    if (!rootedFile.is_open())
        throw std::logic_error("ImportFromFileAndCanonicalize: rootedFile not open!\n");
    count = 0;
    while(std::getline(rootedFile, tempLine))
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)

        GraphContainer temp(n, m, g6String, true, 2); /// NOTE: rooted graphs should be canonical
        temp.SetRootedVertex(0, 0);
        temp.SetRootedVertex(1, 1);
#ifdef DEBUG
        if (temp.GetL()!=l)
            throw std::invalid_argument("ImportFromFileAndCanonicalize requires the correct number of bonds!\n");
#endif
        this->AddCanonicalRootedGraph(temp);
        count++;
    }
    rootedFile.close();
#ifdef DEBUG
    std::cout << "IMPORTED_ROOTED: (" << n << "," << l << "): " << count << "\n";
#endif

}

void CanonicalGraphManager::AddCanonicalGraph(const GraphContainer& g)
{
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    this->CanonicalGraphs[g.GetL()-1].push_back(g);
    this->NbrGraphs++;
}

void CanonicalGraphManager::AddCanonicalRootedGraph(const GraphContainer& g)
{
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    this->CanonicalRootedGraphs[g.GetL()-1].push_back(g);
    this->NbrRootedGraphs++;
}

int CanonicalGraphManager::GetNbrGraphs(int l) const
{
    if (l<1 || l>this->LMax)
        throw std::invalid_argument("GetNbrGraphs requires 1 <= l <= LMax!\n");
    return this->CanonicalGraphs[l-1].size();
}

int CanonicalGraphManager::GetNbrRootedGraphs(int l) const
{
    if (l<1 || l>this->LMax)
        throw std::invalid_argument("GetNbrRootedGraphs requires 1 <= l <= LMax!\n");
    return this->CanonicalRootedGraphs[l-1].size();
}
