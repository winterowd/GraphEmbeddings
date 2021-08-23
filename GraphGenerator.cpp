#include "GraphGenerator.h"

#include <iostream>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

GenericConnectedUndirectedGraphGenerator::GenericConnectedUndirectedGraphGenerator(int argc, char *argv[]) :
    Parameters(argc, argv),
    Graph(this->Parameters.GetN(), this->Parameters.GetLMax())
{}

void GenericConnectedUndirectedGraphGenerator::Generate(bool verbose)
{
    std::vector<std::vector<int>> tab(this->GetN(), std::vector<int>(this->GetLMax(), 0));

    std::vector<int> s(this->Graph.GetLMax(), -1);
    std::vector<int> colK(this->Graph.GetLMax(), -1);

    int count = 1;
    int sk, k, row, col; /// start at ZERO!

    this->AddBond(0,1); /// first graph
    s[0] = 0;
    colK[0] = 1;
    std::cout << "NEW_GRAPH FOUND! number " << count << "\n" << Graph << "\n";

    /// variable initialization
    k = 1;
    s[1] = 0;
    colK[1] =2;

    while (1) /// main loop
    {
        sk = s[k] + 1;
        s[k] = sk;

        col = this->GetColForGenerate(sk); /// column for new bond

        if (sk>=this->Graph.GetNTimesNMinusOneDiv2() || (col-colK[k-1])>1) /// check if we have reached the last hole or graph is disconnected/isolated vertices
        {
            if (verbose)
            {
                if ((col-colK[k-1])>1)
                    std::cout << "DISCONNECTED!\n";
                if (sk>=this->Graph.GetNTimesNMinusOneDiv2())
                    std::cout << "REACHED LAST HOLE!\n";
            }
            k--; /// backtrack
            if (k==0) /// never need to remove the first bond! Largest bit in canonical code! GENERATION COMPLETE!
                break;
            /// reset variables
            sk = s[k];
            row = this->GetRowForGenerate(sk);
            col = colK[k];
            if (verbose)
                std::cout << "Backtrack RemoveBond " << k+1 << " between " << row+1 << " and " << col+1 << "\n";
            this->RemoveBond(row, col); /// remove bond
        }
        else /// connected and have not reached the last hole
        {

            bool backtrackAfterInnerLoop = false; /// do backtracking step after exiting inner while loop?
            bool continueInner = true; /// break out of inner while loop?

            while(continueInner)
            {
                row = this->GetRowForGenerate(sk); /// row for new bond
                if (verbose)
                    std::cout << "AddBond " << k+1 << " between " << row+1 << " and " << col+1 << "\n";
                this->AddBond(row, col); /// add bond

                colK[k] = col; /// update column k

                if (this->IsCanonical(col, verbose)) /// check if new graph is canonical
                {
                    // new graph produced
                    count++;
                    std::cout << "NEW_GRAPH FOUND! number " << count << "\n" << this->Graph << "\n";

                    tab[col][k]++;

                    k++; /// increment k

                    if (k>=this->GetLMax()) /// reached max number of bonds
                    {
                        continueInner = false; /// break out of inner while loop
                        backtrackAfterInnerLoop = true; /// backtrack!
                    }
                    else
                    {
                        sk = s[k-1]+1;

                        if(sk>=this->Graph.GetNTimesNMinusOneDiv2()) /// we have reached the last hole
                        {
                            //std::cout << "reached last hole!\n";
                            continueInner = false; /// exit inner while loop
                            backtrackAfterInnerLoop = true; /// backtrack!
                        }
                        else /// continue inside inner while loop
                        {
                            s[k] = sk;
                            col = this->GetColForGenerate(sk);
                            colK[k] = col;
                        }
                    }
                }
                else /// not canonical; remove bond and break out of inner while loop
                {
                    if (verbose)
                        std::cout << "Not Canonical RemoveBond " << k+1 << " between " << row+1 << " and " << col+1 << "\n";
                    this->RemoveBond(row, col);
                    continueInner = false;
                }

                if (backtrackAfterInnerLoop)
                {
                    if (verbose)
                        std::cout << "backtrack " << k << "\n";
                    k--; /// backtrack

                    // reset variables
                    sk = s[k];

                    row = this->GetRowForGenerate(sk);
                    col = colK[k];

                    if (verbose)
                        std::cout << "Backtrack RemoveBond " << k+1 << " between " << row+1 << " and " << col+1 << "\n";
                    this->RemoveBond(row, col); /// remove bond
                }

            } /// inner while

        } /// else

    } /// outer while
    std::cout << "TOTAL_GRAPHS: " << count << "\n";

    for (int i=2; i<tab.size(); ++i)
        for (int j=0; j<tab[i].size(); ++j)
            std::cout << "TAB: vertices " << i+1 << " links " << j+1 << " " << tab[i][j] << "\n";
}

bool GraphParameters::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                (",n", po::value<unsigned int>(&this->N)->required(), "Maximum vertex order")
                (",l", po::value<unsigned int>(&this->LMax)->required(), "Maximum number of bonds")
                (",d", po::value<bool>(&this->Disconnected)->default_value(false), "Allow disconnected graphs?")
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

GraphParameters::GraphParameters(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}
