#include <iostream>
#include <string>

#include "GraphGeneratorNauty.h"

extern "C" {
int
GENG_MAIN(int argc, char *argv[]);
}

#include <boost/program_options.hpp>
namespace po = boost::program_options;

GraphGeneratorNauty::GraphGeneratorNauty(int argc, char *argv[]) : Parameters(argc, argv)
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

        GENG_MAIN(argv.size()-1,argv.data());

        arguments.pop_back();
        arguments.pop_back();

    }

}

bool GraphGeneratorParametersNauty::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                (",n", po::value<unsigned int>(&this->N)->required(), "Maximum vertex order")
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

GraphGeneratorParametersNauty::GraphGeneratorParametersNauty(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}
