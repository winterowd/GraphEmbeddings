#include "CanonicalGraphManager.h"

/// constructor
/// @param lMax: maximum number of bonds for container
CanonicalGraphManager::CanonicalGraphManager(int lMax) :
    LMax(lMax),
    NbrGraphs(0),
    NbrRootedGraphs(0),
    CanonicalGraphs(this->LMax, std::vector<GraphContainer>()),
    CanonicalRootedGraphs(this->LMax, std::vector<GraphContainer>())
{
    this->GenerateCanonicalGraphs();
}

/// call GetGraphsNauty for all values (n,l) such that we get all graphs with links 1 <= L <= LMax
void CanonicalGraphManager::GenerateCanonicalGraphs()
{
    for (int n=2; n<=this->LMax+1; ++n) ///
    {
        int connUpper = n*(n-1)/2;
        int minBondsFixedN = n-1;
        int maxBondsFixedN = (connUpper < this->LMax) ? connUpper : this->LMax;
#ifdef DEBUG
        std::cout << "DEBUG_GENERATE N: " << n << " L: [" << minBondsFixedN << "," << maxBondsFixedN << "]\n";
#endif
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

    arguments.push_back("Foo"); /// dummy first argument
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
/// @param n: number of vertices
/// @param l: number of links
/// @param filenameUnrooted: filename for unrooted graphs (standard name from GraphGeneratorNauty object)
/// @param filenameRooted: filename for rooted graphs (standard name from GraphGeneratorNauty object)
void CanonicalGraphManager::ImportFromFileAndCanonicalize(int n, int l, std::string filenameUnrooted, std::string filenameRooted)
{
    std::ifstream unrootedFile(filenameUnrooted, std::ifstream::in); /// unrooted file
    if (!unrootedFile.is_open())
        throw std::logic_error("ImportFromFileAndCanonicalize: unrootedFile not open!\n");
    std::string tempLine;
    int m = SETWORDSNEEDED(n);
    int count=0;
    while (std::getline(unrootedFile, tempLine)) /// loop over unrooted graphs (canonicalize before we add)
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)
#ifdef DEBUG
        auto temp = AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(n, g6String);
        if (temp.GetL()!=l || temp.GetN()!=n)
            throw std::invalid_argument("ImportFromFileAndCanonicalize requires the correct number of bonds and vertices!\n");
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
    while(std::getline(rootedFile, tempLine)) /// loop over rooted graphs (already canonical)
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)

        GraphContainer temp(n, m, g6String, true, 2); /// NOTE: rooted graphs should be canonical
        temp.SetRootedVertex(0, 0);
        temp.SetRootedVertex(1, 1);
#ifdef DEBUG
        if (temp.GetL()!=l || temp.GetN()!=n)
            throw std::invalid_argument("ImportFromFileAndCanonicalize requires the correct number of bonds and vertices!\n");
#endif
        this->AddCanonicalRootedGraph(temp);
        count++;
    }
    rootedFile.close();
#ifdef DEBUG
    std::cout << "IMPORTED_ROOTED: (" << n << "," << l << "): " << count << "\n";
#endif

}

/// add a canonical graph to container
/// @param g: canonical unrooted graph container
void CanonicalGraphManager::AddCanonicalGraph(const GraphContainer& g)
{
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    this->CanonicalGraphs[g.GetL()-1].push_back(g);
    this->NbrGraphs++;
}

/// add a canonical rooted graph to container
/// @param g: canonical rooted graph container
void CanonicalGraphManager::AddCanonicalRootedGraph(const GraphContainer& g)
{
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    this->CanonicalRootedGraphs[g.GetL()-1].push_back(g);
    this->NbrRootedGraphs++;
}

/// get the number of graphs with a certain number of bonds
/// @param l: number of bonds
int CanonicalGraphManager::GetNbrGraphs(int l) const
{
    if (l<1 || l>this->LMax)
        throw std::invalid_argument("GetNbrGraphs requires 1 <= l <= LMax!\n");
    return this->CanonicalGraphs[l-1].size();
}

/// get the number of rooted graphs with a certain number of bonds
/// @param l: number of bonds
int CanonicalGraphManager::GetNbrRootedGraphs(int l) const
{
    if (l<1 || l>this->LMax)
        throw std::invalid_argument("GetNbrRootedGraphs requires 1 <= l <= LMax!\n");
    return this->CanonicalRootedGraphs[l-1].size();
}

/// get the 2D index of a given unrooted graph: (L, graphIndex)
/// NOTE: expects container to be already canonical with respect to vertex labels (produced by NAUTY)
/// @param container: GraphContainer containing the graph whose index we want
std::pair<int, int> CanonicalGraphManager::GetGraphIndex(const GraphContainer& container)
{
    if (container.GetL() > this->LMax)
        throw std::invalid_argument("GetGraphIndex requires container's number of bonds to satisfy 1 <= L <= LMax!\n");
    for (int i=0; i<this->GetNbrGraphs(container.GetL()); ++i)
    {
        if (container==this->CanonicalGraphs[container.GetL()-1][i])
            return std::pair<int,int>(container.GetL(),i);
    }
#ifdef DEBUG
    std::cerr << "WARNING! Requested graph not found in list of canonical unrooted graphs!\n";
#endif
    return std::pair<int,int>(-1,-1);
}

/// get the 2D index of a given rooted graph: (L, graphIndex)
/// NOTE: expects container to be already canonical with respect to vertex labels (produced by NAUTY)
/// @param container: GraphContainer containing the graph whose index we want
std::pair<int, int> CanonicalGraphManager::GetRootedGraphIndex(const GraphContainer& container)
{
    if (container.GetL() > this->LMax)
        throw std::invalid_argument("GetRootedGraphIndex requires container's number of bonds to satisfy 1 <= L <= LMax!\n");
    for (int i=0; i<this->GetNbrRootedGraphs(container.GetL()); ++i)
        if (container==this->CanonicalRootedGraphs[container.GetL()-1][i])
            return std::pair<int,int>(container.GetL(),i);
    std::cerr << "WARNING! Requested rooted graph not found in list of canonical rooted graphs!\n";
    return std::pair<int,int>(-1,-1);
}

/// accessor for unrooted graph
/// @param l: number of bonds
/// @param graphIndex: index of the graph in the container for a given number of bonds
GraphContainer CanonicalGraphManager::GetGraph(int l, int graphIndex) const
{
    if (l < 1 || l > this->LMax)
        throw std::invalid_argument("GetGraph expects number of bonds to satisfy 1 <= L <= L_Max!\n");
    if (graphIndex < 0 || graphIndex >= this->CanonicalGraphs[l-1].size())
        throw std::invalid_argument("GetGraph expects graph number to satisfy 0 <= graphNbr < N_graphs at order L!\n");
    return this->CanonicalGraphs[l-1][graphIndex];
}

/// accessor for rooted graph
/// @param l: number of bonds
/// @param graphIndex: index of the graph in the container for a given number of bonds
GraphContainer CanonicalGraphManager::GetRootedGraph(int l, int graphIndex) const
{
    if (l < 1 || l > this->LMax)
        throw std::invalid_argument("GetRootedGraph expects number of bonds to satisfy 1 <= L <= L_Max!\n");
    if (graphIndex < 0 || graphIndex >= this->CanonicalGraphs[l-1].size())
        throw std::invalid_argument("GetRootedGraph expects graph number to satisfy 0 <= graphNbr < N_graphs at order L!\n");
    return this->CanonicalRootedGraphs[l-1][graphIndex];
}
