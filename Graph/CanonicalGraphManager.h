#ifndef CANONICALGRAPHMANAGER_H
#define CANONICALGRAPHMANAGER_H

#include <iostream>
#include <vector>

#include "GraphContainer.h"
#include "GraphGeneratorNauty.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// TODO: takes in L_max and generates all connected graphs from L=1,..,L_Max (canonicalizes in the process)
class CanonicalGraphManager
{
private:
    /**** private variables ****/
    int LMax;

    int NbrGraphs;

    int NbrRootedGraphs;

    std::vector<std::vector<GraphContainer>> CanonicalGraphs;

    std::vector<std::vector<GraphContainer>> CanonicalRootedGraphs;

    /**** private methods ****/
    void GenerateCanonicalGraphs();

    void GetGraphsNauty(int n, int l);

    void ImportFromFileAndCanonicalize(int n, int l, std::string filenameUnrooted, std::string filenameRooted);

    void AddCanonicalGraph(const GraphContainer& g);

    void AddCanonicalRootedGraph(const GraphContainer& g);

    bool FileExists(const std::string& name);

public:
    CanonicalGraphManager(int lMax, bool rooted=false);

    /**** accessors ****/

    int GetTotalNbrGraphs() const { return this->NbrGraphs; }

    int GetTotalNbrRootedGraphs() const { return this->NbrRootedGraphs; }

    int GetNbrGraphs(int l) const;

    int GetNbrRootedGraphs(int l) const;

    GraphContainer GetGraph(int l, int index) const;

    GraphContainer GetRootedGraph(int l, int index) const;

};

inline bool CanonicalGraphManager::FileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

#endif // CANONICALGRAPHMANAGER_H
