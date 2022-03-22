#ifndef MYLAMBDAPOLYNOMIAL_H
#define MYLAMBDAPOLYNOMIAL_H

#include <ginac/ginac.h>

#include "AbstractLattice.h"
#include "AuxiliaryRoutinesForGINAC.h"

/// basic wrapper for polynomial
template<typename T>
class MyLambdaPolynomial
{
private:
    GiNaC::ex Polynomial; /// GiNaC expression

    int MaxManhattanDistance; /// maximum Manhattan distance

public:
    MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance);

    GiNaC::ex GetPolynomial() const { return this->Polynomial; } /// accessor
};

#endif // MYLAMBDAPOLYNOMIAL_H
