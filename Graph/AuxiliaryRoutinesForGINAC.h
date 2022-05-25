#ifndef AUXILIARYROUTINESFORGNAC_H
#define AUXILIARYROUTINESFORGNAC_H

#include <ginac/ginac.h>

#include "AbstractLattice.h"

namespace AuxiliaryRoutinesForGinac {

const GiNaC::symbol& GetSymbol(int l, bool lambda=true);
std::pair<std::string, std::string> MaxInteractionLengthToSymbolName(int l);
std::pair<std::string, std::string> GetFermionCoupling(int index);
int ComputeManhattanDistance(const std::array<int, MaxInteractionLength::NbrInteractions>& coefficients);
GiNaC::ex GetLambdaExpandedPureGaugeRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsDenom, int maxManhattanDistance);
GiNaC::ex GetLambdaExpandedSingleFlavorRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsDenom, int maxManhattanDistance);
GiNaC::ex GetLambdaExpandedRationalFunction(GiNaC::ex num, GiNaC::ex denom, int maxManhattanDistance);
GiNaC::ex GetSingleFlavorH1AndHBar1ExpandedFunction(GiNaC::ex ex, int maxOrderH1, int maxOrderHBar1);
GiNaC::ex GetFullyExpandedSingleFlavorRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsDenom, int maxManhattanDistance, int maxOrderH1, int maxOrderHBar1);
GiNaC::ex GetFullyExpandedSingleFlavorRationalFunction(GiNaC::ex num, GiNaC::ex denom, int maxManhattanDistance, int maxOrderH1, int maxOrderHBar1);
template<typename T>
T CreateRationalPolynomialCoefficient(int num, int denom);

/// create the polynomial in \lambda_i's from the coefficients
/// @param coefficients: vector of coefficients for a given order defined by n_i's
/// @param maxManhattanDistance: maximum Manhattan distance which we use to check if coefficients are valid (debug mode)
template<typename T>
inline GiNaC::ex PureGaugePolynomialFromCoefficients(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance)
{
#ifdef DEBUG
    for (int i=0; i<coefficients.size(); ++i)
        if (AuxiliaryRoutinesForGinac::ComputeManhattanDistance(coefficients[i].second)>maxManhattanDistance)
            std::cerr << "WARNING: AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients encountered a coefficient corresponding to a term larger than the requested Manhattan distance!\n";
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

/// create the polynomial in \lambda_i's and {h_1, \bar{h}_1} from the coefficients
/// @param coefficients: vector of coefficients for a given order defined by n_i's
/// @param maxManhattanDistance: maximum Manhattan distance which we use to check if coefficients are valid (debug mode)
template<typename T>
inline GiNaC::ex SingleFlavorPolynomialFromCoefficients(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficients, int maxManhattanDistance)
{
#ifdef DEBUG
    for (int i=0; i<coefficients.size(); ++i)
        if (AuxiliaryRoutinesForGinac::ComputeManhattanDistance(coefficients[i].second)>maxManhattanDistance)
            std::cerr << "WARNING: AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients encountered a coefficient corresponding to a term larger than the requested Manhattan distance!\n";
#endif
    GiNaC::ex result;
    for (int i=0; i<coefficients.size(); ++i)
    {
        GiNaC::ex tempTerm = coefficients[i].first;
        for (int j=0; j<MaxInteractionLength::NbrInteractions; ++j) /// powers of \lambda_i's
            tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(j), coefficients[i].second[j]);
        /// t^{MD}, MD: Manhattan distance corresponding to given product of \lambda_i's
        std::array<int, MaxInteractionLength::NbrInteractions> lambdaCoefficients;
        for (int j=0; j<MaxInteractionLength::NbrInteractions; ++j)
            lambdaCoefficients[j] = coefficients[i].second[j];
        tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(MaxInteractionLength::NbrInteractions), AuxiliaryRoutinesForGinac::ComputeManhattanDistance(lambdaCoefficients));
        /// powers of h_1, \bar{h}_1
        /// TODO: multiple flavors
        tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(0, false), coefficients[i].second[MaxInteractionLength::NbrInteractions]);
        tempTerm *= GiNaC::pow(AuxiliaryRoutinesForGinac::GetSymbol(1, false), coefficients[i].second[MaxInteractionLength::NbrInteractions+1]);
        result += tempTerm;
    }
    return result;
}

}

#endif // AUXILIARYROUTINESFORGNAC_H
