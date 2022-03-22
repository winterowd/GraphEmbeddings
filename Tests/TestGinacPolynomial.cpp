#include "AuxiliaryRoutinesForGINAC.h"
#include "MyLambdaPolynomial.h"

//GiNaC::ex PolynomialFromCoefficients(const std::vector<std::pair<double, std::array<int, MaxInteractionLength::NbrInteractions>>>& coefficients, int maxManhattanDistance);

int main()
{
    std::vector<std::pair<double, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients1;
    std::vector<std::pair<GiNaC::numeric, std::array<int, MaxInteractionLength::NbrInteractions>>> coefficients2;
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder1{1,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder2{3,0,0,0};
    std::array<int, MaxInteractionLength::NbrInteractions> tempOrder3{1,1,0,0};

    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1., tempOrder1});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{3., tempOrder2});
    coefficients1.push_back(std::pair<double,std::array<int, MaxInteractionLength::NbrInteractions>>{1./6., tempOrder3});

    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1), tempOrder1});
    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(3), tempOrder2});
    coefficients2.push_back(std::pair<GiNaC::numeric,std::array<int, MaxInteractionLength::NbrInteractions>>{GiNaC::numeric(1,6), tempOrder3});

    auto result1 = AuxiliaryRoutinesForGinac::PolynomialFromCoefficients<double>(coefficients1, 6);
    auto result2 = AuxiliaryRoutinesForGinac::PolynomialFromCoefficients<GiNaC::numeric>(coefficients2, 6);

    std::cout << "RESULT1: " << result1 << "\n";
    std::cout << "RESULT2: " << result2 << "\n";

    auto result3 = MyLambdaPolynomial<double>(coefficients1, 6);

    return 0;
}
