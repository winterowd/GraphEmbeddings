#include "MyLambdaPolynomial.h"

/// constructor
template<typename T>
MyLambdaPolynomial<T>::MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance) :
    Polynomial(AuxiliaryRoutinesForGinac::PolynomialFromCoefficients(coefficients, maxManhattanDistance)),
    MaxManhattanDistance(maxManhattanDistance)
{

}

/// only need to template on double or GiNaC::numeric (should be an integer or a fraction!)
template class MyLambdaPolynomial<double>;
template class MyLambdaPolynomial<GiNaC::numeric>;
