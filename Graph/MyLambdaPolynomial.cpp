#include "MyLambdaPolynomial.h"

/// constructor (pure gauge)
template<typename T>
MyLambdaPolynomial<T>::MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance) :
    Polynomial(AuxiliaryRoutinesForGinac::PureGaugePolynomialFromCoefficients(coefficients, maxManhattanDistance)),
    MaxManhattanDistance(maxManhattanDistance)
{

}

/// constructor (single flavor of static quarks)
template<typename T>
MyLambdaPolynomial<T>::MyLambdaPolynomial(const std::vector<std::pair<T, std::array<int, MaxInteractionLength::NbrInteractions+2>>>& coefficients, int maxManhattanDistance) :
    Polynomial(AuxiliaryRoutinesForGinac::SingleFlavorPolynomialFromCoefficients(coefficients, maxManhattanDistance)),
    MaxManhattanDistance(maxManhattanDistance)
{

}

template<typename T>
MyLambdaPolynomial<T>& MyLambdaPolynomial<T>::operator+=(const MyLambdaPolynomial<T>& rhs)
{
    if (this->MaxManhattanDistance!=rhs.GetMaxManhattanDistance())
        throw std::invalid_argument("MyLambdaPolynomial::operator+= requires rhs to have the same MaxManhattanDistance!\n");
    this->Polynomial += rhs.GetPolynomial();
    return *this;
}

template<typename T>
MyLambdaPolynomial<T> operator+(const MyLambdaPolynomial<T>& p1, const MyLambdaPolynomial<T>& p2)
{
    if (p1.GetMaxManhattanDistance()!=p2.GetMaxManhattanDistance())
        throw std::invalid_argument("MyLambdaPolynomial::operator+ requires p1 and p2 to have the same MaxManhattanDistance!\n");
    return MyLambdaPolynomial<T>(p1.GetPolynomial()+p2.GetPolynomial(), p1.GetMaxManhattanDistance());
}

template<typename T>
MyLambdaPolynomial<T> operator+(const MyLambdaPolynomial<T>& p, const GiNaC::ex& value)
{
    return MyLambdaPolynomial<T>(p.GetPolynomial()+value, p.GetMaxManhattanDistance());
}

template<typename T>
MyLambdaPolynomial<T> operator+(const GiNaC::ex& value, const MyLambdaPolynomial<T>& p)
{
    return MyLambdaPolynomial<T>(p.GetPolynomial()+value, p.GetMaxManhattanDistance());
}

/// only need to template on double or GiNaC::numeric (should be an integer or a fraction!)
template class MyLambdaPolynomial<double>;
template class MyLambdaPolynomial<GiNaC::numeric>;
