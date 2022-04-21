#ifndef AUXILIARYROUTINESFORGNAC_H
#define AUXILIARYROUTINESFORGNAC_H

#include <ginac/ginac.h>

#include "AbstractLattice.h"

namespace AuxiliaryRoutinesForGinac {

const GiNaC::symbol& GetSymbol(int l);
std::pair<std::string, std::string> MaxInteractionLengthToSymbolName(int l);
int ComputeManhattanDistance(const std::array<int, MaxInteractionLength::NbrInteractions>& coefficients);
GiNaC::ex GetExpandedRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsDenom, int maxManhattanDistance);
GiNaC::ex GetExpandedRationalFunction(GiNaC::ex num, GiNaC::ex denom, int maxManhattanDistance);
template<typename T>
T CreateRationalPolynomialCoefficient(int num, int denom);

/// create the polynomial in \lambda_i's from the coefficients
/// @param coefficients: vector of coefficients for a given order defined by n_i's
/// @param
/// @param maxManhattanDistance: maximum Manhattan distance which we use to check if coefficients are valid (debug mode)
template<typename T>
inline GiNaC::ex PolynomialFromCoefficients(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance)
{
#ifdef DEBUG
    for (int i=0; i<coefficients.size(); ++i)
        if (AuxiliaryRoutinesForGinac::ComputeManhattanDistance(coefficients[i].second)>maxManhattanDistance)
            std::cerr << "WARNING: AuxiliaryRoutinesForGinac::PolynomialFromCoefficients encountered a coefficient corresponding to a term larger than the requested Manhattan distance!\n";
#endif
    GiNaC::ex result;
    for (int i=0; i<coefficients.size(); ++i)
    {
        GiNaC::ex tempTerm = coefficients[i].first;
        for (int j=0; j<MaxInteractionLength::NbrInteractions; ++j) /// powers of \lambda_i's
            tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(j), coefficients[i].second[j]);
        /// t^{MD}, MD: Manhattan distance corresponding to given product of \lambda_i's
        tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(MaxInteractionLength::NbrInteractions), AuxiliaryRoutinesForGinac::ComputeManhattanDistance(coefficients[i].second));
        result += tempTerm;
    }
    return result;
}

}

#endif // AUXILIARYROUTINESFORGNAC_H
