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

    bool PureGauge; /// pure gauge?

public:
    /// first constructor (pure gauge)
    MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance);

    /// second constructor (single flavor of static quarks)
    MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficients, int maxManhattanDistance);

    /// third constructor (generic)
    MyLambdaPolynomial(const GiNaC::ex& polynomial, int maxManhattanDistance) : Polynomial(polynomial), MaxManhattanDistance(maxManhattanDistance) {}

    GiNaC::ex GetPolynomial() const { return this->Polynomial; } /// accessor

    int GetMaxManhattanDistance() const { return this->MaxManhattanDistance; } /// accessor

    bool IsPureGauge() const { return this->PureGauge; }

    /// overloaded operator
    MyLambdaPolynomial<T>& operator+=(const MyLambdaPolynomial<T>& rhs);

};

/// overloaded operators
template<typename T>
MyLambdaPolynomial<T> operator+(const MyLambdaPolynomial<T>& p1, const MyLambdaPolynomial<T>& p2);
template<typename T>
MyLambdaPolynomial<T> operator+(const MyLambdaPolynomial<T>& p, const GiNaC::ex& value);
template<typename T>
MyLambdaPolynomial<T> operator+(const GiNaC::ex& value, const MyLambdaPolynomial<T>& p);

#endif // MYLAMBDAPOLYNOMIAL_H
