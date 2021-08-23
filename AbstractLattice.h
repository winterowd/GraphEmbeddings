#ifndef ABSTRACTLATTICE_H
#define ABSTRACTLATTICE_H

#include <vector>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

enum MaxInteractionLength
{
    NearestNeighbor,
    NextNearestNeighbor,
    NbrInteractions /// change this to number of lengths available and use in other files
};

inline std::istream& operator>>(std::istream &in, MaxInteractionLength& l)
{
    std::string token;
    in >> token;
    std::cout << "TOKEN: " << token << "\n";
    if (token == "NN")
    {
        l = MaxInteractionLength::NearestNeighbor;
    }
    else if (token == "NNN")
    {
        l = MaxInteractionLength::NextNearestNeighbor;
    }
    else
    {
        throw po::validation_error(po::validation_error::invalid_option_value, "max interaciton length");
    }

    return in;
}

class AbstractLattice {

private:

    std::string Name;

    MaxInteractionLength Length;

public:

    AbstractLattice(std::string name, MaxInteractionLength length=NearestNeighbor) : Name(name), Length(length) {}

    virtual ~AbstractLattice(){};

    virtual int GetSiteIndex(const std::vector<unsigned int>& indices) const  = 0;

    virtual void GetSiteCoordinates(int index, std::vector<unsigned int>& indices) = 0;

    virtual int GetNearestNeighbor(int siteIndex, int nnIndex) = 0;

    virtual int GetNextNearestNeighbor(int siteIndex, int nnnIndex) = 0;

    virtual unsigned int GetN() const = 0;

    virtual unsigned int GetNbrNN() const = 0;

    virtual unsigned int GetNbrNNN() const = 0;

    virtual unsigned int GetDim() const = 0;

    virtual bool AreNN(unsigned int index1, unsigned int index2) = 0;

    virtual bool AreNNN(unsigned int index1, unsigned int index2) = 0;

    std::string GetName() const { return this->Name; }

};

#endif // ABSTRACTLATTICE_H
