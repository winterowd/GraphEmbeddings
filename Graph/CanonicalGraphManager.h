#ifndef CANONICALGRAPHMANAGER_H
#define CANONICALGRAPHMANAGER_H

#include <iostream>
#include <vector>

#include "GraphContainer.h"
#include "GraphGeneratorNauty.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// takes in L_max and generates all connected unrooted AND rooted graphs from L=1,..,L_Max (canonicalizes in the process)
class CanonicalGraphManager
{
private:
    /**** private variables ****/
    int LMax; /// maximum number of links (>1)

    int NbrGraphs; /// total number of unrooted graphs

    std::vector<int> NbrRootedGraphs; /// total number of rooted graphs (one and two roots)

    std::vector<std::vector<GraphContainer>> CanonicalGraphs; /// container for unrooted graphs

    std::vector<std::vector<GraphContainer>> CanonicalOneRootedGraphs; /// container for one-rooted graphs

    std::vector<std::vector<GraphContainer>> CanonicalTwoRootedGraphs; /// container for two-rooted graphs

    /**** private methods ****/
    void GenerateCanonicalGraphs();

    void GetGraphsNauty(int n, int l);

    void ImportUnrootedFromFileAndCanonicalize(int n, int l, std::string filenameUnrooted);

    void ImportRootedFromFileAndCanonicalize(int n, int l, std::string filenameRooted, int nbrRoots);

    void AddCanonicalGraph(const GraphContainer& g);

    void AddCanonicalRootedGraph(const GraphContainer& g, int nbrRoots);

    bool FileExists(const std::string& name);

    GraphContainer GetOneRootedGraph(int l, int graphIndex) const;

    GraphContainer GetTwoRootedGraph(int l, int graphIndex) const;

    std::pair<int, int> GetUnrootedGraphIndex(const GraphContainer& container);

    std::pair<int, int> GetRootedGraphIndex(const GraphContainer& container, int nbrRoots);

public:
    CanonicalGraphManager(int lMax);

    /**** public accessors ****/

    int GetLMax() const { return this->LMax; }

    int GetTotalNbrGraphs() const { return this->NbrGraphs; }

    int GetTotalNbrRootedGraphs(int nbrRoots) const;

    int GetNbrGraphs(int l) const;

    int GetNbrRootedGraphs(int l, int nbrRoots) const;

    std::pair<int, int> GetGraphIndex(const GraphContainer& container);

    GraphContainer GetGraph(int l, int graphIndex) const;

    GraphContainer GetRootedGraph(int l, int graphIndex, int nbrRoots) const;

};

inline bool CanonicalGraphManager::FileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

#endif // CANONICALGRAPHMANAGER_H
