

#include "GraphGeneratorNauty.h"

extern "C" {
int
GENG_MAIN(int argc, char *argv[]);
}

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/// constructor which receives command line arguments from user
GraphGeneratorNauty::GraphGeneratorNauty(int argc, char *argv[]) :
    Parameters(argc, argv),
    fp(NULL)
{}

/// Set up geng argument list.  The 0-th argument is the command name. There must be a NULL at the end.
void GraphGeneratorNauty::Generate()
{

    std::vector<std::string> arguments;

    arguments.push_back("geng");
    if (!this->Parameters.AllowDisconnected())
        arguments.push_back("-c");

    for (int i=2; i<=this->Parameters.GetN(); ++i)
    {
        arguments.push_back(std::to_string(i));

        std::string filename = "graphs_g6_";
        if (!this->Parameters.AllowDisconnected())
            filename += "connected_";
        filename += "order_"+std::to_string(i)+".dat";

        arguments.push_back(filename);

        std::vector<char*> argv;
        for (const auto& arg: arguments)
            argv.push_back((char*)arg.data());
        argv.push_back(nullptr);

        /// call routine which reads in graphs in filename and then generates rooted graphs
        if (this->Parameters.GenerateTwoRooted())
        {
            /// compare speeds of routines to generate rooted graphs (debugging)
            auto start = std::chrono::high_resolution_clock::now();
            this->GenerateTwoRootedFixedOrder(i, filename, this->Parameters.IsVerbose(), true);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> s_int = end-start;
            std::cout << "FIXED_ORDER: " << i << " " << s_int.count() << "s\n";
            start = std::chrono::high_resolution_clock::now();
            this->GenerateTwoRootedFixedOrderIterative(i, filename, this->Parameters.IsVerbose(), true);
            end = std::chrono::high_resolution_clock::now();
            s_int = end-start;
            std::cout << "FIXED_ORDER_ITERATIVE: " << i << " " << s_int.count() << "s\n";
        }
        else /// otherwise generate connected graphs
            GENG_MAIN(argv.size()-1,argv.data());

        arguments.pop_back();
        arguments.pop_back();

    }

}

/// test routine for relabeling
void GraphGeneratorNauty::TestRelabeling(int n, std::string inputFilename)
{
    this->N = n;

    this->GetAllPossiblePairsForRootedVertices();

    /// read in a random graph (should match order n)
    int m = SETWORDSNEEDED(this->N);

    FILE *fp = fopen(inputFilename.c_str(), "r");

    if (fp!=NULL)
        std::cout << "TestRelabeling: Successfully opened " << inputFilename << std::endl;
    else
        throw std::invalid_argument("TestRelabeling: Error opening "+inputFilename);

    DYNALLSTAT(graph, g, g_sz);
    DYNALLOC2(graph, g, g_sz, n, m, "malloc");

    graph *gtemp = readg(fp,g,m,&m,&n);

    GraphContainer refGraph(this->N, m, g);

    std::vector<int> newLabels(this->N, -1);
    for (int i=0; i<this->RootedVertexNumbers.size(); ++i)
    {
        this->ProduceNewLabelingGivenRootedVertices(this->RootedVertexNumbers[i], newLabels, true);
        GraphContainer tempGraph(refGraph);
        tempGraph.RelabelVertices(newLabels);
    }

    DYNFREE(g,g_sz);
    fclose(fp);
}

/// set the vertex colors based on an array of rooted vertices
/// vertex label rootedVertices[i] given color i, where i=0,1,...,k-1 (assume all elements of rootedVertices to be unique!)
/// all other vertices given color k (size of rootedVertices)
/// c: C-style array of N elements for colors (assume it is allocated elsewhere!)
/// rootedVertices: vector containing labels of rooted vertices
void GraphGeneratorNauty::SetVertexColors(int *c, const std::vector<int>& rootedVertices, bool verbose)
{
    if (rootedVertices.size()>this->N)
        throw std::invalid_argument("SetVertexColors requires number of rooted vertices to be less than or equal to N!\n");

    if (std::find_if(rootedVertices.begin(), rootedVertices.end(), [this](const int& x) { return x>=this->N; } ) != rootedVertices.end())
        throw std::invalid_argument("SetVertexColors requires each rooted vertex label x to satisfy 0 <= x < N!\n");

    if (verbose)
        std::cout << "SetVertexColors: " << rootedVertices.size() << " rooted vertices!\n";

    for (int i=0; i<this->N; ++i)
        c[i] = rootedVertices.size(); /// all vertices given color k

    for (int i=0; i<rootedVertices.size(); ++i) /// rooted vertices given colors 0,1,...,k-1
    {
        c[rootedVertices[i]] = i;
        if (verbose)
            std::cout << "SetVertexColors: vertex " << rootedVertices[i] << " assigned color " << i << "\n";
    }
}

/// sets up the nauty data structures lab and ptn according to our vertex coloring convention (N vertices with 3 colors, rooted vertices being colors 0 and 1 and ALWAYS assigned labels 0 and 1!!!!!)
/// c: C-style array of N integers specifying the color of each vertex i.e. c_j gives the color of vertex labeled j (should be 0, 1, or 2)
/// lab: C-style array of N integers specifying vertex labels
/// ptn: C-style array of N integers specifying and partition of the vertex labels into colors
void GraphGeneratorNauty::SetColoredPartition(int* c, int* lab, int* ptn)
{
    /* loop over all colors to fill it */
    int cur_index = 0;
    for (int i = 0; i<3; i++) /// loop over colors
    {
        for (int j = 0; j<this->N; j++) /// loop over vertex labels
        {
            if (c[j] == i)
            {
                lab[cur_index] = j;
                ptn[cur_index] = 1;
                cur_index++;
            }
        }

        if (cur_index > 0) /// partition ends with a zero
        {
            ptn[cur_index - 1] = 0;
        }
    }
}

/// generate two-rooted graphs from list of connected graphs of order N
/// read in graphs which are assumed to be generated by geng and enumerate the two-rooted graphs
/// generate all one-rooted graphs using the help of densenauty
/// this is done by choosing each of N vertices, canonicalizing via densenauty and adding to a list, making sure that no duplicates are added
/// generate all two-rooted graphs using the previously generated list of one-rooted graphs
/// this is done by choosing each of remaining N-1 vertices (vertex 0 always labeled color 0 by convention), canonicalizing via densenauty and adding to list (BEHIND 1-rooted graphs)
/// n: number of vertices
/// inputFilename: name of the file containing the graphs in g6 format
/// verbose: flag for verbosity
/// outputSorted: sort output for comparison
void GraphGeneratorNauty::GenerateTwoRootedFixedOrderIterative(int n, std::string inputFilename, bool verbose, bool outputSorted)
{
    int countComparison = 0; /// for debugging number of comparisons

    /// set up auxiliary variables for NAUTY
    this->N = n; /// set graph size
    this->MWords = SETWORDSNEEDED(this->N); /// set number of words

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canong

    nauty_check(WORDSIZE, this->N, this->MWords, NAUTYVERSIONID); /// check if everything is ok

    this->fp = fopen(inputFilename.c_str(), "r"); /// open file containing graphs of a fixed order from which we will generate two-rooted graphs
    if (this->fp!=NULL)
        std::cout << "GenerateTwoRootedFixedOrder: Successfully opened " << inputFilename << std::endl;
    else
        throw std::invalid_argument("GenerateTwoRootedFixedOrder: Error opening "+inputFilename);

    //// output file stream for rooted graphs
    std::string outputFilename = "test_iterative_graphs_g6_rooted_";
    if (!this->Parameters.AllowDisconnected())
        outputFilename += "connected_";
    outputFilename += "order_"+std::to_string(this->N)+".dat";
    FILE *fpg6 = fopen(outputFilename.c_str(), "w");

    FILE *fpg6sorted;
    if (outputSorted)
    {
        std::string outputSortedFilename = "sorted_test_iterative_graphs_g6_rooted_";
        if (!this->Parameters.AllowDisconnected())
            outputSortedFilename += "connected_";
        outputSortedFilename += "order_"+std::to_string(this->N)+".dat";
        fpg6sorted = fopen(outputSortedFilename.c_str(), "w");
    }

    /// declare nauty data structures
    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg1, cg1_sz); /// declare canonical graphs
    DYNALLSTAT(graph, cg2, cg2_sz); /// first carries the canonical relabeling after the first rooted vertex chosen and the second after both rooted vertices are chosen
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    /// allocate nauty data structures
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(graph, cg1, cg1_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(graph, cg2, cg2_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, lab, lab_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, this->N, this->MWords, "malloc");

    int *c = (int*)malloc(this->N * sizeof(int)); /// alloc C-style arrray for colors of vertices

    GraphContainer refContainer(this->N, this->MWords, true, 2); /// container for graphs
    int count = 0; /// initialize counter
    while (1) /// read graphs in file (ASSUME THEY ARE ALL OF THE SAME ORDER N!)
    {
        graph *g1 = this->GetNextGraph(g);
        if (g1 == NULL)
            break;
        count++; /// increment counter

        std::vector<GraphContainer> rootedGraphList; /// holds the list of rooted graphs produced from a single connected graph (one line of file produced by geng)

        refContainer.SetGraphFromDenseNauty(g); /// setup container from nauty dense format
        int countComparisonIndividual = 0; /// number of comparisons to create all two-rooted graphs for each connected, simple graph produced by geng
        for (int i=0; i<this->N; ++i) /// loop over first vertex
        {
            GraphContainer tempGraph(refContainer); /// copy constructor

            tempGraph.SetRootedVertex(0,i); /// set first rooted vertex

            /// set colors
            std::vector<int> tempRootedList(1,i);
            this->SetVertexColors(c, tempRootedList);

            this->SetColoredPartition(c, lab, ptn); /// set coloring (will be written over in call to densenauty)

            /// call densenauty
            densenauty(g, lab, ptn, orbits, &options, &stats, this->MWords, this->N, cg1);

            tempGraph.ColoredCanonicalRelabeling(lab, i); /// relabel according to canonical labeling returned from densenauty (first rooted vertex only)
#ifdef DEBUG
            GraphContainer testGraphFromDenseCG(this->N, this->MWords, cg1, true, 2);
            testGraphFromDenseCG.SetRootedVertex(0,0);
            if (testGraphFromDenseCG!=tempGraph)
                throw std::invalid_argument("In GenerateTwoRootedFixedOrderIterative relabeled graph does not equal graph created from canonical dense NAUTY graph afer first rooted vertex placed!\n");
#endif
            /// add tempGraph to rootedGraphList if not already in the list
            auto tempIt = std::find(rootedGraphList.begin(), rootedGraphList.end(), tempGraph);
            if (tempIt == rootedGraphList.end())
            {
                rootedGraphList.push_back(tempGraph);
                if (verbose)
                {
                    countComparison += rootedGraphList.size();
                    countComparisonIndividual += rootedGraphList.size();
                    std::cout << "ADDING ROOTED GRAPH TO LIST (FIRST ROOTED VERTEX)!\n";
                }
            }
            else
            {
                if (verbose)
                {
                    countComparison += std::distance(rootedGraphList.begin(), tempIt)+1;
                    countComparisonIndividual += std::distance(rootedGraphList.begin(), tempIt)+1;
                    std::cout << "DUPLICATE! NOT ADDING TO LIST (FIRST ROOTED VERTEX)!\n";
                }
            }
        } /// loop over first rooted vertices

        auto sizeAfterFirstSweep = rootedGraphList.size(); /// list extends past the one-rooted graphs
        for (int i=0; i<sizeAfterFirstSweep; ++i) /// loop over one-rooted graphs in list
        {
            rootedGraphList[i].GetDenseNautyFromGraph(cg1); /// get dense nauty form of current graph
            for (int j=1; j<this->N; ++j) /// loop over vertices
            {
                GraphContainer tempGraph(rootedGraphList[i]); /// copy constructor (using graph already in list)

                tempGraph.SetRootedVertex(1,j); /// set second rooted vertex

                /// set colors
                std::vector<int> tempRootedList{0, j};
                this->SetVertexColors(c, tempRootedList);

                this->SetColoredPartition(c, lab, ptn); /// set coloring (will be written over in call to densenauty)

                /// call densenauty
                densenauty(cg1, lab, ptn, orbits, &options, &stats, this->MWords, this->N, cg2);

                tempGraph.ColoredCanonicalRelabeling(lab, 0, j); /// relabel according to canonical labeling returned from densenauty
#ifdef DEBUG
                GraphContainer testGraphFromDenseCG(this->N, this->MWords, cg2, true, 2);
                testGraphFromDenseCG.SetRootedVertex(0,0);
                testGraphFromDenseCG.SetRootedVertex(1,1);
                if (testGraphFromDenseCG!=tempGraph)
                    throw std::invalid_argument("In GenerateTwoRootedFixedOrderIterative relabeled graph does not equal graph created from canonical dense NAUTY graph afer second rooted vertex placed!\n");
#endif
                /// add tempGraph to rootedGraphList if not already in the list (start AFTER the elements added on the first sweep)
                auto tempIt = std::find(rootedGraphList.begin()+sizeAfterFirstSweep, rootedGraphList.end(), tempGraph);
                if (tempIt == rootedGraphList.end())
                {
                    rootedGraphList.push_back(tempGraph);
                    if (verbose)
                    {
                        countComparison += rootedGraphList.size()-sizeAfterFirstSweep;
                        countComparisonIndividual += rootedGraphList.size()-sizeAfterFirstSweep;
                        std::cout << "ADDING ROOTED GRAPH TO LIST (SECOND ROOTED VERTEX)!\n";
                    }
                    /// output cg to file
                    char *s = ntog6(cg2,this->MWords,this->N);
                    fprintf(fpg6, "%.*s ", static_cast<int>(strlen(s)-1), s); /// g6 string to file without newline character
                    if (!outputSorted)
                        writegroupsize(fpg6,stats.grpsize1,stats.grpsize2); /// write group size (symmetry factor) to file
                    fprintf(fpg6,"\n");
                }
                else
                    if (verbose)
                    {
                        countComparison += std::distance(rootedGraphList.begin()+sizeAfterFirstSweep, tempIt)+1;
                        countComparisonIndividual += std::distance(rootedGraphList.begin()+sizeAfterFirstSweep, tempIt)+1;
                        std::cout << "DUPLICATE! NOT ADDING TO LIST (SECOND ROOTED VERTEX)!\n";
                    }
            } /// loop over vertices

        } /// loop over graphs in list after first sweep

        if (verbose)
        {
            std::cout << "At order " << this->N << " graph number " << count << " produced " << rootedGraphList.size() << " rooted graphs!\n";
            char *s = ntog6(g,this->MWords,this->N);
            std::cout << "TOTAL_COUNTS_INDIVIDUAL_ITERATIVE: " << this->N << " " << count << " " << countComparisonIndividual << " " << s << "\n";
        }

        if (outputSorted) /// want a sorted list to compare output of two procedures
        {
            std::sort(rootedGraphList.begin()+sizeAfterFirstSweep, rootedGraphList.end()); /// sort list of rooted graphs produced from a single simple connected graph
            for (int i=sizeAfterFirstSweep; i<rootedGraphList.size(); ++i)
            {
                rootedGraphList[i].GetDenseNautyFromGraph(cg1); /// get dense nauty format
                char *s = ntog6(cg1,this->MWords,this->N); /// convert to g6 string
                fprintf(fpg6sorted, "%.*s \n", static_cast<int>(strlen(s)-1), s); /// output so that diff'ing will work!
            }
        }

    } /// while

    if (verbose)
        std::cout << "TOTAL_COUNTS_ITERATIVE at order " << this->N << ": " << countComparison << "\n";

    /// free nauty data structures
    DYNFREE(g, g_sz);
    DYNFREE(cg1, cg1_sz);
    DYNFREE(cg2, cg2_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    free(c); /// free vertex colors

    fclose(fpg6); /// close output file
    if (outputSorted)
        fclose(fpg6sorted);
    fclose(fp); /// close input file

}

/// generate two-rooted graphs from list of connected graphs of order N
/// read in graphs which are assumed to be generated by geng and enumerate the two-rooted graphs
/// generate all pairs of vertices which will be colors 0 and 1 (rooted graphs). Relabel the graph so that the two rooted vertices will have labels 0 and 1, respectively
/// call densenauty to get canonical ordering of vertices of color 2 (rooted vertices map to themselves)
/// compare adjacency matrices of previously generated rooted graphs coming from a single unrooted graph to see if we have generated a new rooted graph
/// n: number of vertices
/// inputFilename: name of the file containing the graphs in g6 format
/// verbose: flag for verbosity
/// outputSorted: sort output for comparison
void GraphGeneratorNauty::GenerateTwoRootedFixedOrder(int n, std::string inputFilename, bool verbose, bool outputSorted)
{
    int countComparison = 0; /// for debugging number of comparisons

    /// set up auxiliary variables for NAUTY
    this->N = n; /// set graph size
    this->MWords = SETWORDSNEEDED(this->N); /// set number of words

    static DEFAULTOPTIONS_GRAPH(options); /// options
    options.defaultptn = false; /// color the vertices
    options.getcanon = true; /// get canong

    nauty_check(WORDSIZE, this->N, this->MWords, NAUTYVERSIONID); /// check if everything is ok

    this->fp = fopen(inputFilename.c_str(), "r"); /// open file containing graphs of a fixed order from which we will generate two-rooted graphs
    if (this->fp!=NULL)
        std::cout << "GenerateTwoRootedFixedOrder: Successfully opened " << inputFilename << std::endl;
    else
        throw std::invalid_argument("GenerateTwoRootedFixedOrder: Error opening "+inputFilename);

    //// output file stream for rooted graphs
    std::string outputFilename = "graphs_g6_rooted_";
    if (!this->Parameters.AllowDisconnected())
        outputFilename += "connected_";
    outputFilename += "order_"+std::to_string(this->N)+".dat";
    FILE *fpg6 = fopen(outputFilename.c_str(), "w");

    FILE *fpg6sorted;
    if (outputSorted)
    {
        std::string outputSortedFilename = "sorted_graphs_g6_rooted_";
        if (!this->Parameters.AllowDisconnected())
            outputSortedFilename += "connected_";
        outputSortedFilename += "order_"+std::to_string(this->N)+".dat";
        fpg6sorted = fopen(outputSortedFilename.c_str(), "w");
    }

    /// generate all pairs of vertices given order N
    this->GetAllPossiblePairsForRootedVertices();

    /// declare nauty data structures
    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLSTAT(graph, cg, cg_sz); /// declare canonical graph
    DYNALLSTAT(int, lab, lab_sz); /// label
    DYNALLSTAT(int, ptn, ptn_sz); /// partition for coloring
    DYNALLSTAT(int, orbits, orbits_sz); /// orbits when calling densenauty
    statsblk stats; /// status

    /// allocate nauty data structures
    DYNALLOC2(graph, g, g_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(graph, cg, cg_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, lab, lab_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, ptn, ptn_sz, this->N, this->MWords, "malloc");
    DYNALLOC2(int, orbits, orbits_sz, this->N, this->MWords, "malloc");

    int *c = (int*)malloc(this->N * sizeof(int)); /// alloc C-style arrray for colors of vertices

    GraphContainer refContainer(this->N, this->MWords); /// container for graphs
    int count = 0; /// initialize counter
    while (1) /// read graphs in file (ASSUME THEY ARE ALL OF THE SAME ORDER N!)
    {
        graph *g1 = this->GetNextGraph(g);
        if (g1 == NULL)
            break;
        count++; /// increment counter

        std::vector<GraphContainer> rootedGraphList; /// holds the list of rooted graphs produced from a single connected graph (one line of file produced by geng)

        refContainer.SetGraphFromDenseNauty(g); /// setup container from nauty dense format
        int countComparisonIndividual = 0; /// number of comparisons to create all two-rooted graphs for each connected, simple graph produced by geng
        for (int i=0; i<this->RootedVertexNumbers.size(); ++i) /// loop over pairs of rooted vertices
        {
            GraphContainer tempGraph(refContainer); /// copy constructor

            /// set colors from RootedVertexNumbers[i]
            this->SetVertexColors(c, this->RootedVertexNumbers[i]);

            this->SetColoredPartition(c, lab, ptn); /// set coloring (will be written over in call to densenauty)

            /// call densenauty
            densenauty(g, lab, ptn, orbits, &options, &stats, this->MWords, this->N, cg);

            tempGraph.ColoredCanonicalRelabeling(lab, this->RootedVertexNumbers[i][0], this->RootedVertexNumbers[i][1]); /// relabel according to canonical labeling returned from densenauty
#ifdef DEBUG
            GraphContainer testGraphFromDenseCG(this->N, this->MWords, cg);
            if (testGraphFromDenseCG!=tempGraph)
                throw std::invalid_argument("In GenerateTwoRootedFixedOrder relabeled graph does not equal graph created from canonical dense NAUTY graph!\n");
#endif
            /// add tempGraph to rootedGraphList if not already in the list
            auto tempIt = std::find(rootedGraphList.begin(), rootedGraphList.end(), tempGraph);
            if (tempIt == rootedGraphList.end())
            {
                rootedGraphList.push_back(tempGraph);
                if (verbose)
                {
                    countComparison += rootedGraphList.size();
                    countComparisonIndividual += rootedGraphList.size();
                    std::cout << "ADDING ROOTED GRAPH TO LIST!\n";
                }
                /// output cg to file
                char *s = ntog6(cg,this->MWords,this->N);
                fprintf(fpg6, "%.*s ", static_cast<int>(strlen(s)-1), s); /// g6 string to file without newline character
                if (!outputSorted)
                    writegroupsize(fpg6,stats.grpsize1,stats.grpsize2); /// write group size (symmetry factor) to file
                fprintf(fpg6,"\n");
            }
            else
                if (verbose)
                {
                    countComparison += std::distance(rootedGraphList.begin(), tempIt)+1;
                    countComparisonIndividual += std::distance(rootedGraphList.begin(), tempIt)+1;
                    std::cout << "DUPLICATE! NOT ADDING TO LIST!\n";
                }

        } /// loop over pairs

        if (verbose)
        {
            std::cout << "At order " << this->N << " graph number " << count << " produced " << rootedGraphList.size() << " rooted graphs!\n";
            char *s = ntog6(g,this->MWords,this->N);
            std::cout << "TOTAL_COUNTS_INDIVIDUAL_NONITERATIVE: " << this->N << " " << count << " " << countComparisonIndividual << " " << s << "\n";
        }

        if (outputSorted) /// want a sorted list to compare output of two procedures
        {
            std::sort(rootedGraphList.begin(), rootedGraphList.end()); /// sort list of rooted graphs produced from a single simple connected graph
            for (int i=0; i<rootedGraphList.size(); ++i)
            {
                rootedGraphList[i].GetDenseNautyFromGraph(cg); /// get dense nauty format
                char *s = ntog6(cg,this->MWords,this->N); /// convert to g6 string
                fprintf(fpg6sorted, "%.*s \n", static_cast<int>(strlen(s)-1), s); /// output so that diff'ing will work!
            }
        }

    } /// while

    if (verbose)
        std::cout << "TOTAL_COUNTS_NONITERATIVE at order " << this->N << ": " << countComparison << "\n";

    /// free nauty data structures
    DYNFREE(g, g_sz);
    DYNFREE(cg, cg_sz);
    DYNFREE(lab,lab_sz);
    DYNFREE(ptn,ptn_sz);
    DYNFREE(orbits,orbits_sz);

    free(c); /// free vertex colors

    fclose(fpg6); /// close output file
    if (outputSorted)
        fclose(fpg6sorted);
    fclose(fp); /// close input file

}

/// given a selected pair vertices, produce a new relabeling
/// rooted vertices should be mapped to 0 and 1 respectively (start from zero!)
void GraphGeneratorNauty::ProduceNewLabelingGivenRootedVertices(const std::vector<int>& rooted, std::vector<int>& newLabeling, bool verbose)
{
    if (rooted.size()!=2)
        throw std::invalid_argument("ProduceNewLabelingGivenRootedVertices expects rooted to have two elements!\n");

    if (newLabeling.size()!=this->N)
        throw std::invalid_argument("ProduceNewLabelingGivenRootedVertices expects newLabeling to have N elements!\n");

    for (int i=0; i<newLabeling.size(); ++i)
        newLabeling[i] = i;

    if (rooted[0]==1 && rooted[1]==0) /// simply swap labels of two rooted vertices?
    {
        newLabeling[0] = 1;
        newLabeling[1] = 0;
    }
    else
    {
        if (rooted[0]!=0) /// does the first rooted vertex have the correct label? (0) Swap!
        {
            newLabeling[0] = rooted[0]; /// vertex 0 (which is in place) receives new label
            newLabeling[rooted[0]] = 0; /// relabel vertex of color 0 as number 0
        }

        if (rooted[1]!=1) /// does the second rooted vertex have the correct label? (1) Swap!
        {
            // find where in array "1" resides and swap places
            auto element = std::find(newLabeling.begin(), newLabeling.end(), 1);
            *element = newLabeling[rooted[1]]; /// vertex 1 receives new label
            newLabeling[rooted[1]] = 1; /// relabel vertex of color 1 as number 1
        }
    }

    if (verbose)
    {
        std::cout << "RELABELING: " << rooted[0] << " " << rooted[1] << "\n";
        for (int i=0; i<newLabeling.size(); ++i)
            std::cout << i << " maps to " << newLabeling[i] << "\n";
    }

}

/// got all pairs of vertices to be fixed for two-point function (labels start at zero!)
/// assumes we have N vertices
void GraphGeneratorNauty::GetAllPossiblePairsForRootedVertices()
{
    std::vector<int> tmp;
    std::vector<int> vertices(this->N);
    for (int i = 0; i<vertices.size(); ++i)
        vertices[i] = i;
    this->RootedVertexNumbers.clear(); /// clear list of fixed vertices
    this->GenerateUniqueCombinationsWithNoDuplicates(tmp, vertices, 2);
}

/// get all unique combinations of vertices without any duplicate entries using depth first search
void GraphGeneratorNauty::GenerateUniqueCombinationsWithNoDuplicates(std::vector<int>& tmp, const std::vector<int>& vertices, int k, bool verbose)
{
    if (k==0)
    {
        if (verbose)
        {
            std::cout << "ADDING_PAIR:";
            for (auto x:tmp)
                std::cout << " " << x;
            std::cout << std::endl;
        }
        this->RootedVertexNumbers.push_back(tmp);
        return;
    }

    for (int i=0; i<vertices.size(); ++i)
    {
        if (std::find(tmp.begin(), tmp.end(), vertices[i])==tmp.end())
        {
            tmp.push_back(vertices[i]);
            this->GenerateUniqueCombinationsWithNoDuplicates(tmp, vertices, k-1);
            tmp.pop_back();
        }
    }
}

/// wrapper for readg
graph* GraphGeneratorNauty::GetNextGraph(graph *g)
{
    return readg(this->fp,g,0,&this->MWords,&this->N);
}

/// process commmand line arguments using BOOST
bool GraphGeneratorParametersNauty::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                (",n", po::value<unsigned int>(&this->N)->required(), "Maximum vertex order")
                (",d", po::value<bool>(&this->Disconnected)->default_value(false), "Allow disconnected graphs?")
                (",r", po::value<bool>(&this->TwoRooted)->default_value(false), "Generate two-rooted graphs?")
                (",v", po::value<bool>(&this->Verbose)->default_value(false), "Verbose?")
                ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help"))
        {
            std::cout << desc << "\n";
            return false;
        }

        po::notify(vm);
    }
    catch(std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    catch (...)
    {
        std::cerr << "Unknown error!" << "\n";
        return false;
    }
    return true;
}

/// constructor for user parameters wrapper class
GraphGeneratorParametersNauty::GraphGeneratorParametersNauty(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}
