#include "CanonicalGraphManager.h"

/// constructor
/// @param lMax: maximum number of bonds for container
CanonicalGraphManager::CanonicalGraphManager(int lMax) :
    LMax(lMax),
    NbrGraphs(0),
    NbrRootedGraphs(2 /* one- and two-rooted graphs only */, 0),
    CanonicalGraphs(this->LMax, std::vector<GraphContainer>()),
    CanonicalOneRootedGraphs(this->LMax, std::vector<GraphContainer>()),
    CanonicalTwoRootedGraphs(this->LMax, std::vector<GraphContainer>())
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

    /// import rooted and unrooted graphs from files
    this->ImportUnrootedFromFileAndCanonicalize(n, l, generator.GetOutputFilenameFixedOrder());
    this->ImportRootedFromFileAndCanonicalize(n, l, generator.GetOutputFilenameFixedOrder(true, 1), 1);
    this->ImportRootedFromFileAndCanonicalize(n, l, generator.GetOutputFilenameFixedOrder(true, 2), 2);

}

/// read in g6 from file: all graphs (unrooted) have n vertices and l bonds
/// @param n: number of vertices
/// @param l: number of links
/// @param filenameUnrooted: filename for unrooted graphs (standard name from GraphGeneratorNauty object)
void CanonicalGraphManager::ImportUnrootedFromFileAndCanonicalize(int n, int l, std::string filenameUnrooted)
{
    std::ifstream unrootedFile(filenameUnrooted, std::ifstream::in); /// unrooted file
    if (!unrootedFile.is_open())
        throw std::logic_error("ImportFromFileAndCanonicalize: unrootedFile not open!\n");
    std::string tempLine;
    int count=0;
    while (std::getline(unrootedFile, tempLine)) /// loop over unrooted graphs (canonicalize before we add)
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)

        this->AddCanonicalGraph(AuxiliaryRoutinesForNauty::GetCanonicalGraphNauty(n, g6String));

        count++;
    }
    unrootedFile.close();
#ifdef DEBUG
    std::cout << "IMPORTED_UNROOTED: (" << n << "," << l << "): " << count << "\n";
#endif
}

/// read in g6 from file: all graphs (rooted) have n vertices and l bonds
/// @param n: number of vertices
/// @param l: number of links
/// @param filenameRooted: filename for rooted graphs (standard name from GraphGeneratorNauty object)
/// @param nbrRoots: number of roots (should match with filenameRooted)
void CanonicalGraphManager::ImportRootedFromFileAndCanonicalize(int n, int l, std::string filenameRooted, int nbrRoots)
{
    std::ifstream rootedFile(filenameRooted, std::ifstream::in); /// unrooted file
    if (!rootedFile.is_open())
        throw std::logic_error("ImportFromFileAndCanonicalize: rootedFile not open!\n");
    std::string tempLine;
    int m = SETWORDSNEEDED(n);
    int count = 0;
    while(std::getline(rootedFile, tempLine)) /// loop over rooted graphs (already canonical)
    {
        std::string g6String;
        std::stringstream ss(tempLine);

        ss >> g6String; /// first token is g6 string (regardless of rooted or unrooted)

        GraphContainer temp(n, m, g6String, nbrRoots); /// NOTE: rooted graphs should be canonical
        /// if canonical, vertices 0 (and 1) should have colors 0 (and 1), where we count vertices from 0,...,N-1
        for (int i=0; i<nbrRoots; ++i)
            temp.SetRootedVertex(i, i);

        /// get the symmetry factor
        int symmFactor;
        ss >> symmFactor;
        temp.SetSymmFactor(symmFactor);
#ifdef DEBUG
        if (temp.GetL()!=l || temp.GetN()!=n)
            throw std::invalid_argument("ImportFromFileAndCanonicalize requires the correct number of bonds and vertices!\n");
        std::vector<int> rootedVertices(nbrRoots);
        for (int i=0; i<nbrRoots; ++i)
            rootedVertices[i] = i;
        auto tempCanonical = AuxiliaryRoutinesForNauty::GetCanonicalColoredGraphNauty(n, g6String, rootedVertices);
        if (tempCanonical!=temp)
            std::cerr << "WARNING: ImportFromFileAndCanonicalize found rooted graph of order " << nbrRoots << " which is not canonical!\n";
#endif
        this->AddCanonicalRootedGraph(temp, nbrRoots);
        count++;
    }
    rootedFile.close();
#ifdef DEBUG
    std::cout << "IMPORTED_ROOTED_" << nbrRoots << ": (" << n << "," << l << "): " << count << "\n";
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
void CanonicalGraphManager::AddCanonicalRootedGraph(const GraphContainer& g, int nbrRoots)
{
    if (nbrRoots != 1 && nbrRoots !=2)
        throw std::invalid_argument("AddCanonicalRootedGraph expects nbrRoots to be 1 or 2!\n");
    if (g.GetL() < 1 || g.GetL() > this->LMax)
        throw std::invalid_argument("AddCanonicalRootedGraph expects g to have number of bonds satisfying 1 <= L <= L_Max!\n");
    if (nbrRoots==1)
        this->CanonicalOneRootedGraphs[g.GetL()-1].push_back(g);
    if (nbrRoots==2)
        this->CanonicalTwoRootedGraphs[g.GetL()-1].push_back(g);
    //std::cout << "DEBUG_AddCanonicalRootedGraph: " << nbrRoots << " vs " << this->NbrRootedGraphs.size() << " vs " << this->N
    this->NbrRootedGraphs[nbrRoots-1]++;
}

/// get the total number of rooted graphs for a given number of rooted vertices
/// @param nbrRoots: number of rooted vertices (1 or 2)
int CanonicalGraphManager::GetTotalNbrRootedGraphs(int nbrRoots) const
{
    if (nbrRoots!=1 && nbrRoots!=2)
        throw std::invalid_argument("GetTotalNbrRootedGraphs requires 1 <= l <= LMax!\n");
    return this->NbrRootedGraphs[nbrRoots-1];
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
int CanonicalGraphManager::GetNbrRootedGraphs(int l, int nbrRoots) const
{
    if (nbrRoots != 1 && nbrRoots !=2)
        throw std::invalid_argument("GetNbrRootedGraphs expects nbrRoots to be 1 or 2!\n");
    if (l<1 || l>this->LMax)
        throw std::invalid_argument("GetNbrRootedGraphs requires 1 <= l <= LMax!\n");
    if (nbrRoots==1)
        return this->CanonicalOneRootedGraphs[l-1].size();
    else
        return this->CanonicalTwoRootedGraphs[l-1].size();
}

/// public accessor for graphs (rooted AND unrooted)
/// NOTE: expects container to be already canonical with respect to vertex labels (produced by NAUTY)
/// @param container: the graph whose index we want
std::pair<int, int> CanonicalGraphManager::GetGraphIndex(const GraphContainer& container)
{
    if (container.GetL() > this->LMax)
        throw std::invalid_argument("GetGraphIndex requires container's number of bonds to satisfy 1 <= L <= LMax!\n");

    if (container.GetNbrRooted()==0)
        return this->GetUnrootedGraphIndex(container);
    else
        return this->GetRootedGraphIndex(container, container.GetNbrRooted());
}

/// get the 2D index of a given unrooted graph: (L, graphIndex)
/// NOTE: expects container to be already canonical with respect to vertex labels (produced by NAUTY)
/// @param container: GraphContainer containing the graph whose index we want
std::pair<int, int> CanonicalGraphManager::GetUnrootedGraphIndex(const GraphContainer& container)
{
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
std::pair<int, int> CanonicalGraphManager::GetRootedGraphIndex(const GraphContainer& container, int nbrRoots)
{
    if (nbrRoots != 1 && nbrRoots !=2)
        throw std::invalid_argument("GetRootedGraphIndex expects nbrRoots to be 1 or 2!\n");
    if (container.GetL() > this->LMax)
        throw std::invalid_argument("GetRootedGraphIndex requires container's number of bonds to satisfy 1 <= L <= LMax!\n");
    for (int i=0; i<this->GetNbrRootedGraphs(container.GetL(), nbrRoots); ++i)
    {
        if (nbrRoots==1)
        {
            if (container==this->CanonicalOneRootedGraphs[container.GetL()-1][i])
                return std::pair<int,int>(container.GetL(),i);
        }
        else
        {
            if (container==this->CanonicalTwoRootedGraphs[container.GetL()-1][i])
                return std::pair<int,int>(container.GetL(),i);
        }
    }
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
/// @param nbrRoots: number of rooted vertices
GraphContainer CanonicalGraphManager::GetRootedGraph(int l, int graphIndex, int nbrRoots) const
{
    if (nbrRoots != 1 && nbrRoots !=2)
        throw std::invalid_argument("GetRootedGraph expects nbrRoots to be 1 or 2!\n");
    if (nbrRoots==1)
        return this->GetOneRootedGraph(l, graphIndex);
    else
        return this->GetTwoRootedGraph(l, graphIndex);
}

/// accessor for one-rooted graphs
/// @param l: number of bonds
/// @param graphIndex: index of the graph in the container for a given number of bonds
GraphContainer CanonicalGraphManager::GetOneRootedGraph(int l, int graphIndex) const
{
    if (l < 1 || l > this->LMax)
        throw std::invalid_argument("GetOneRootedGraph expects number of bonds to satisfy 1 <= L <= L_Max!\n");
    if (graphIndex < 0 || graphIndex >= this->CanonicalOneRootedGraphs[l-1].size())
        throw std::invalid_argument("GetOneRootedGraph expects graph number to satisfy 0 <= graphNbr < N_graphs at order L!\n");
    return this->CanonicalOneRootedGraphs[l-1][graphIndex];
}

/// accessor for two-rooted graphs
/// @param l: number of bonds
/// @param graphIndex: index of the graph in the container for a given number of bonds
GraphContainer CanonicalGraphManager::GetTwoRootedGraph(int l, int graphIndex) const
{
    if (l < 1 || l > this->LMax)
        throw std::invalid_argument("GetTwoRootedGraph expects number of bonds to satisfy 1 <= L <= L_Max!\n");
    if (graphIndex < 0 || graphIndex >= this->CanonicalTwoRootedGraphs[l-1].size())
        throw std::invalid_argument("GetTwoRootedGraph expects graph number to satisfy 0 <= graphNbr < N_graphs at order L!\n");
    return this->CanonicalTwoRootedGraphs[l-1][graphIndex];
}
