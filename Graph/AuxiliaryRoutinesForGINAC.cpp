#include "AuxiliaryRoutinesForGINAC.h"

/// compute the Manhattan distance from n_i's (powers of \lambda_i's)
/// ManhattanDistance = n_1 + 2 \times (n_2 + n_3) + 3 \times n_4
/// @param coefficients: set of n_i's
int AuxiliaryRoutinesForGinac::ComputeManhattanDistance(const std::array<int, MaxInteractionLength::NbrInteractions>& coefficients)
{
    return coefficients[0]+2*(coefficients[1]+coefficients[2])+3*coefficients[3];
}

/// TODO: multiple flavors
/// factory for the GiNaC symbol objects \lambda_i and {h_1, \bar{h}_1}
/// this is so perturbative expressions can be combined (see GiNaC manual)
/// @param l: length (0, 1, ..., NbrInteractions-1)
const GiNaC::symbol& AuxiliaryRoutinesForGinac::GetSymbol(int l, bool lambda)
{
    std::pair<std::string, std::string> resultPair;
    if (lambda)
        resultPair = AuxiliaryRoutinesForGinac::MaxInteractionLengthToSymbolName(l);
    else
        resultPair = AuxiliaryRoutinesForGinac::GetFermionCoupling(l);
    //std::cout << "resultPair: " << resultPair.first << " " << resultPair.second << "\n";
    static std::map<std::string, GiNaC::symbol> directory;
    auto it = directory.find(resultPair.first);
    if (it != directory.end())
        return it->second;
    else
    {
        //auto temp = GiNaC::symbol(resultPair.first, resultPair.second);
        //std::cout << "NEW_SYMBOL: " << temp << "\n";
        //auto result = directory.insert(std::make_pair(resultPair.first, GiNaC::symbol(resultPair.first, resultPair.second))).first->second;
        //return result;
        return directory.insert(std::make_pair(resultPair.first, GiNaC::symbol(resultPair.first, resultPair.second))).first->second;
    }
}

/// name and latex name from a length
/// @param l: length (0, 1, ..., NbrInteractions-1)
std::pair<std::string, std::string> AuxiliaryRoutinesForGinac::MaxInteractionLengthToSymbolName(int l)
{
    if (l<0 || l>MaxInteractionLength::NbrInteractions)
        throw std::invalid_argument("MaxInteractionLengthToSymbolName requires a valid l!\n");

    std::pair<std::string, std::string> result;
    if (l==0)
        result = {"lambda_1", "\\lambda_1"};
    else if (l==1)
        result = {"lambda_2", "\\lambda_2"};
    else if (l==2)
        result = {"lambda_3", "\\lambda_3"};
    else if (l==3)
        result = {"lambda_4", "\\lambda_4"};
    else if (l==4) /// MaxInteractionLength::NbrInteractions should be 4!
        result = {"t", "\t"}; /// dummy variable for series expansions

    return result;
}

/// TODO: multiple flavors
/// get the appropriate fermion variable (single flavor)
/// @param index for fermion coupling
std::pair<std::string, std::string> AuxiliaryRoutinesForGinac::GetFermionCoupling(int index)
{
    if (index<0 || index>1)
        throw std::invalid_argument("AuxiliaryRoutinesForGinac::GetFermionCoupling expects index to be 0 or 1!\n");
    std::pair<std::string, std::string> result;
    if (index==0)
        result = {"h_1", "\h_1"};
    else
        result = {"hBar_1", "\\bar{h}_1"};
    return result;
}

/// takes the vector of coefficients, powers of \lambda_i and calls the function below
/// creates the polynomials using our routine PureGaugePolynomialFromCoefficients
/// @param coefficientsNum: coefficients for the numerator
/// @param coefficientsDenom: coefficients for the denominator
/// @param maxManhattanDistance: maximum Manhattan distance
GiNaC::ex AuxiliaryRoutinesForGinac::GetLambdaExpandedPureGaugeRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficientsDenom, int maxManhattanDistance)
{
    return AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients(coefficientsNum, /*true,*/ maxManhattanDistance), AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients(coefficientsDenom, /*true,*/ maxManhattanDistance), maxManhattanDistance);
}

/// function which takes the ratio of two polynomials in variables {\lambda_i, i=0,1,..} and t and expands to a given Manhattan distance (given by power of t)
/// \lim_{t \to 1} f(t)/g(t), where the ratio is expanded about t=0 to a given power
/// @param num: polynomial in the numerator
/// @param denom: polynomial in the denominator
/// @param maxManhattanDistance: maximum Manhattan distance
GiNaC::ex AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(GiNaC::ex num, GiNaC::ex denom, int maxManhattanDistance)
{
    return (num/denom).series(AuxiliaryRoutinesForGinac::GetSymbol(MaxInteractionLength::NbrInteractions)==0, maxManhattanDistance+1).subs(AuxiliaryRoutinesForGinac::GetSymbol(MaxInteractionLength::NbrInteractions) == 1).expand();
}

/// takes the vector of coefficients, powers of \lambda_i and {h_1, \bar{h}_1} and calls the function below
/// creates the polynomials using our routine PureGaugePolynomialFromCoefficients
/// @param coefficientsNum: coefficients for the numerator
/// @param coefficientsDenom: coefficients for the denominator
/// @param maxManhattanDistance: maximum Manhattan distance
GiNaC::ex AuxiliaryRoutinesForGinac::GetLambdaExpandedSingleFlavorRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsDenom, int maxManhattanDistance)
{
    return AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(AuxiliaryRoutinesForGinac::SingleFlavorPolynomialFromCoefficients(coefficientsNum, maxManhattanDistance), AuxiliaryRoutinesForGinac::SingleFlavorPolynomialFromCoefficients(coefficientsDenom, maxManhattanDistance), maxManhattanDistance);
}

/// takes a polynomial and performs an expansion in both h_1 and \hbar{h}_1
/// @param ex: GiNaC expression which should be a polynomial in \{ \lambda_i \} and { h_1, \hbar{h}_1 \}
/// @param maxOrderH1: maximum expansion order for h_1
/// @param maxOrderHBar1: maximum expansion order for \hbar{h}_1
GiNaC::ex AuxiliaryRoutinesForGinac::GetSingleFlavorH1AndHBar1ExpandedFunction(GiNaC::ex ex, int maxOrderH1, int maxOrderHBar1)
{
    auto expandInH1 = GiNaC::series_to_poly(ex.series(AuxiliaryRoutinesForGinac::GetSymbol(0, false)==0, maxOrderH1+1));
    return expandInH1.series(AuxiliaryRoutinesForGinac::GetSymbol(1, false)==0, maxOrderHBar1+1);
}

/// perform a full expansion in  \{ \lambda_i \} h_1, and \hbar{h}_1 (single flavor)
/// @param num: numerator
/// @param denom: denominator
/// @param maxManhattanDistance: maximum Manhattan distance for expansion in \lambda_i
/// @param maxOrderH1: maximum order in h_1
/// @param maxOrderHBar1: maximum order in \hbar{h}_1
GiNaC::ex AuxiliaryRoutinesForGinac::GetFullyExpandedSingleFlavorRationalFunction(GiNaC::ex num, GiNaC::ex denom, int maxManhattanDistance, int maxOrderH1, int maxOrderHBar1)
{
    auto temp = AuxiliaryRoutinesForGinac::GetLambdaExpandedRationalFunction(num, denom, maxManhattanDistance);
    return AuxiliaryRoutinesForGinac::GetSingleFlavorH1AndHBar1ExpandedFunction(temp, maxOrderH1, maxOrderHBar1);
}

/// peform a full expansion in \{ \lambda_i \} h_1, and \hbar{h}_1 (single flavor)
/// @param coefficientsNum: coefficients for the numerator
/// @param coefficientsDenom: coefficients for the denominator
/// @param maxManhattanDistance: maximum Manhattan distance for expansion in \lambda_i
/// @param maxOrderH1: maximum order in h_1
/// @param maxOrderHBar1: maximum order in \hbar{h}_1
GiNaC::ex AuxiliaryRoutinesForGinac::GetFullyExpandedSingleFlavorRationalFunction(const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsNum, const std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficientsDenom, int maxManhattanDistance, int maxOrderH1, int maxOrderHBar1)
{
    auto temp = AuxiliaryRoutinesForGinac::GetLambdaExpandedSingleFlavorRationalFunction(coefficientsNum, coefficientsDenom, maxManhattanDistance);
    return AuxiliaryRoutinesForGinac::GetSingleFlavorH1AndHBar1ExpandedFunction(temp, maxOrderH1, maxOrderHBar1);
}

template<>
double AuxiliaryRoutinesForGinac::CreateRationalPolynomialCoefficient<double>(int num, int denom)
{
    if (denom==0)
        throw std::invalid_argument("AuxiliaryRoutinesForGinac::CreateRationalPolynomialCoefficient<double>: denom is zero!\n");
    return num/denom;
}

template<>
GiNaC::numeric AuxiliaryRoutinesForGinac::CreateRationalPolynomialCoefficient<GiNaC::numeric>(int num, int denom)
{
    if (denom==0)
        throw std::invalid_argument("AuxiliaryRoutinesForGinac::CreateRationalPolynomialCoefficient<GiNaC::numeric>: denom is zero!\n");
    return GiNaC::numeric(num,denom);
}
