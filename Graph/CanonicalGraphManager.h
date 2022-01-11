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

    bool Rooted;

    int TotalGraphs;

    std::vector<std::vector<GraphContainer>> CanonicalGraphs;

    /**** private methods ****/
    void GenerateCanonicalGraphs();

    void RunGraphGeneratorNauty(int n, int l);

    void ImportFromFileAndCanonicalize(int n, int l, std::string filename);

    void AddCanonicalGraph(const GraphContainer& g);

public:
    CanonicalGraphManager(int lMax, bool rooted=false);

};

#endif // CANONICALGRAPHMANAGER_H
