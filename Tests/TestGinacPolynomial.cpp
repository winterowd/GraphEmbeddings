#include "AuxiliaryRoutinesForGINAC.h"
#include "MyLambdaPolynomial.h"
#include "TestTemplatePureGauge.h"

template<typename T>
T SingleSiteIntegral(int n, int m);

/// binomial coefficient for doubles (recursive)
double TestBinomialCoefficient(int n, int k)
{
#ifdef DEBUG
    if (k > n)
        throw std::invalid_argument("BinomialCoefficient requires k <= n!\n");
    if (k < 0 || n <0)
        throw std::invalid_argument("BinomialCoefficient requires both n and k to be greater than zero!\n");
#endif
    if (k==0 || k==n)
        return 1;
    return (TestBinomialCoefficient(n-1, k-1) + TestBinomialCoefficient(n-1, k));
}

/// factorial for double (recursive)
/// @param n: size
double TestFactorial(int n)
{
#ifdef DEBUG
    if (n<0)
        throw std::invalid_argument("Factorial requires a non-negative argument!\n");
#endif
    return (n==0) || (n==1) ? 1 : n*TestFactorial(n-1);
}

template<>
double SingleSiteIntegral<double>(int n, int m)
{
    int nMinusM = n-m;
    if (nMinusM%3!=0)
        return 0;
    int nMinusMDiv3 = nMinusM/3;
    int jMin = nMinusM<0 ? 0 : nMinusMDiv3;
    int jMax = floor(n/3);
    double result = 0;
    for (int j=jMin; j<=jMax; ++j)
    {
        double bin1 = TestBinomialCoefficient(3*(n-j-nMinusMDiv3+1), n-3*j);
        double bin2 = TestBinomialCoefficient(2*j-nMinusMDiv3,j);
        double nFactorial = TestFactorial(n);
        double mFactorial = TestFactorial(m);
        double denom1 = TestFactorial(n-j-nMinusMDiv3+1);
        double denom2 = TestFactorial(n-j-nMinusMDiv3+2);
        double denom3 = TestFactorial(2*j-nMinusMDiv3);
        result += 2*(nFactorial/denom1)*(mFactorial/denom2)*(bin1*bin2/denom3);
    }
    return result;
}

template<>
GiNaC::numeric SingleSiteIntegral<GiNaC::numeric>(int n, int m)
{
    GiNaC::numeric nG = n;
    GiNaC::numeric mG = m;
    int nMinusM = n-m;
    if (nMinusM%3!=0)
        return GiNaC::numeric(0);
    GiNaC::numeric nMinusMDiv3(nMinusM,3);
    int jMin = nMinusM<0 ? 0 : nMinusMDiv3.to_int();
    int jMax = floor(n/3);
    GiNaC::numeric result = 0;
    for (int j=jMin; j<=jMax; ++j)
    {
        GiNaC::numeric bin1 = GiNaC::binomial(3*(nG-j-nMinusMDiv3+1), nG-3*j);
        GiNaC::numeric bin2 = GiNaC::binomial(2*j-nMinusMDiv3,GiNaC::numeric(j));
        GiNaC::numeric nFactorial = GiNaC::factorial(nG);
        GiNaC::numeric mFactorial = GiNaC::factorial(mG);
        GiNaC::numeric denom1 = GiNaC::factorial(n-j-nMinusMDiv3+1);
        GiNaC::numeric denom2 = GiNaC::factorial(n-j-nMinusMDiv3+2);
        GiNaC::numeric denom3 = GiNaC::factorial(2*j-nMinusMDiv3);
        result += 2*(nFactorial/denom1)*(mFactorial/denom2)*(bin1*bin2/denom3);
    }
    if (!result.is_rational())
        throw std::logic_error("PureGaugeWeight<GiNaC::numeric>::SingleSiteIntegral should give a rational number!\n");
    return result;
}

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

    /*** SOME_TESTS_FOR_GINAC::NUMERIC ***/
    int n = 4;
    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz); /// declare graph
    DYNALLOC2(graph, g, g_sz, n, m, "malloc"); /// allocate graph

    std::string g6String = "Cr"; /// "square" graph
    char *tempg6 = new char[g6String.length()+1];

    std::strcpy(tempg6, g6String.c_str());
    stringtograph(tempg6, g, m); /// g6 string to densenauty

    delete[] tempg6;

    GraphContainer WeightContainer(n, m, g); /// container from densenauty
    std::vector<ExternalPolyakovLoop> rootedVertices{ExternalPolyakovLoop{1,true}, ExternalPolyakovLoop{2,false}};

    TestTemplatePureGauge<double> someTest(WeightContainer, rootedVertices);
    TestTemplatePureGauge<GiNaC::numeric> anotherTest(WeightContainer, rootedVertices);
    PureGaugeWeight yetAnotherTest(WeightContainer, rootedVertices);

    std::cout << someTest.Weight() << " " << anotherTest.Weight() << " " << yetAnotherTest.Weight() << "\n";

    return 0;
}
