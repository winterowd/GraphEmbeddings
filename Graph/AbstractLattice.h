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
    ThirdNearestNeighbor,
    FourthNearestNeighbor,
    NbrInteractions
};

inline std::string StringFromMaxInteractionLength(MaxInteractionLength l)
{
    std::string result;
    if (l==MaxInteractionLength::NearestNeighbor)
        result = "NN";
    else if (l==MaxInteractionLength::NextNearestNeighbor)
        result = "NNN";
    else if (l==MaxInteractionLength::ThirdNearestNeighbor)
        result = "3NN";
    else if (l==MaxInteractionLength::FourthNearestNeighbor)
        result = "4NN";
    else
        throw std::invalid_argument("StringFromMaxInteractionLength requires a valid l!\n");
    return result;
}

inline std::istream& operator>>(std::istream &in, MaxInteractionLength& l)
{
    std::string token;
    in >> token;
    if (token == "NN")
    {
        l = MaxInteractionLength::NearestNeighbor;
    }
    else if (token == "NNN")
    {
        l = MaxInteractionLength::NextNearestNeighbor;
    }
    else if (token == "3NN")
    {
        l = MaxInteractionLength::ThirdNearestNeighbor;
    }
    else if (token == "4NN")
    {
        l = MaxInteractionLength::FourthNearestNeighbor;
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

    virtual int GetThirdNearestNeighbor(int siteIndex, int n3Index) = 0;

    virtual int GetFourthNearestNeighbor(int siteIndex, int n4Index) = 0;

    virtual unsigned int GetN() const = 0;

    virtual unsigned int GetNbrNN() const = 0;

    virtual unsigned int GetNbrNNN() const = 0;

    virtual unsigned int GetNbrThirdNN() const = 0; /// third nearest-neighbor

    virtual unsigned int GetNbrFourthNN() const = 0; /// fourth nearest-neighbor

    virtual unsigned int GetDim() const = 0;

    virtual bool AreNN(unsigned int index1, unsigned int index2) = 0;

    virtual bool AreNNN(unsigned int index1, unsigned int index2) = 0;

    virtual bool AreThirdNN(unsigned int index1, unsigned int index2) = 0;

    virtual bool AreFourthNN(unsigned int index1, unsigned int index2) = 0;

    virtual int GetManhattanDistance(int neighborDepth) = 0;

    std::string GetName() const { return this->Name; }

    MaxInteractionLength GetMaxInteractionLength() const { return this->Length; }

};

#endif // ABSTRACTLATTICE_H
