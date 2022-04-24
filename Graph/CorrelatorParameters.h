#ifndef CORRELATORPARAMETERS_H
#define CORRELATORPARAMETERS_H

#include "AbstractLattice.h"
#include "MiscellaneousRoutines.h"

class CorrelatorParameters /// correlators parameters (get from user using Boost program options)
{
private:
    int MaxManhattanDistance; /// maximum Manhattan distance

    MaxInteractionLength MaxEmbeddingLength; /// length of longest bond for embedding

    MaxInteractionLength CorrelatorLength; /// correlator length

    bool StaticQuarks; /// include static quarks?

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    CorrelatorParameters(int argc, char *argv[]);

    int GetMaxManhattanDistance() const { return this->MaxManhattanDistance; }

    MaxInteractionLength GetMaxEmbeddingLength() const { return this->MaxEmbeddingLength; }

    MaxInteractionLength GetCorrelatorLength() const { return this->CorrelatorLength; }

    bool IncludeStaticQuarks() const { return this->StaticQuarks; }

};

/// process user's command-line arguments
inline bool CorrelatorParameters::ProcessCommandLine(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "Produce help message")
                ("order,n", po::value<int>(&this->MaxManhattanDistance)->required(), "Maximum Manhattan distance")
                (",d", po::value<MaxInteractionLength>(&this->MaxEmbeddingLength)->default_value(MaxInteractionLength::NearestNeighbor), "MaxEmbeddingLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                (",m", po::value<MaxInteractionLength>(&this->CorrelatorLength)->default_value(MaxInteractionLength::NearestNeighbor), "CorrelatorLength: NN NNN 3N or 4N (longer distances not yet supported!)")
                (",s", po::value<bool>(&this->StaticQuarks)->default_value(false), "Static quarks?")
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

/// constructor
inline CorrelatorParameters::CorrelatorParameters(int argc, char *argv[])
{
    if (this->ProcessCommandLine(argc, argv))
        std::cout << "Successfully parsed command line arguments!" << "\n";
    else
        std::exit(1);
}

#endif // CORRELATORPARAMETERS_H
