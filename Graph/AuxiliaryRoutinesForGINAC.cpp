#include "AuxiliaryRoutinesForGINAC.h"

/// compute the Manhattan distance from n_i's (powers of \lambda_i's)
/// ManhattanDistance = n_1 + 2 \times (n_2 + n_3) + 3 \times n_4
/// @param coefficients: set of n_i's
int AuxiliaryRoutinesForGinac::ComputeManhattanDistance(const std::array<int, MaxInteractionLength::NbrInteractions>& coefficients)
{
    return coefficients[0]+2*(coefficients[1]+coefficients[2])+3*coefficients[3];
}

/// factory for the GiNaC symbol objects \lambda_i
/// this is so perturbative expressions can be combined (see GiNaC manual)
/// @param l: length (0, 1, ..., NbrInteractions-1)
const GiNaC::symbol& AuxiliaryRoutinesForGinac::GetSymbol(int l)
{
    auto resultPair = AuxiliaryRoutinesForGinac::MaxInteractionLengthToSymbolName(l);
    static std::map<std::string, GiNaC::symbol> directory;
    auto it = directory.find(resultPair.first);
    if (it != directory.end())
        return it->second;
    else
        return directory.insert(std::make_pair(resultPair.first, GiNaC::symbol(resultPair.first, resultPair.second))).first->second;
}

/// name and latex name from a length
/// @param l: length (0, 1, ..., NbrInteractions-1)
std::pair<std::string, std::string> AuxiliaryRoutinesForGinac::MaxInteractionLengthToSymbolName(int l)
{
    std::pair<std::string, std::string> result;
    if (l==0)
        result = {"lambda_1", "\\lambda_1"};
    else if (l==1)
        result = {"lambda_2", "\\lambda_2"};
    else if (l==2)
        result = {"lambda_3", "\\lambda_3"};
    else if (l==3)
        result = {"lambda_4", "\\lambda_4"};
    else
        throw std::invalid_argument("MaxInteractionLengthToSymbolName requires a valid l!\n");
    return result;
}
