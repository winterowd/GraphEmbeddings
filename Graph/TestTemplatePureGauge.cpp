#include "TestTemplatePureGauge.h"

/// constructor
/// @arg container: pointer to GraphContainer object describing the graph we are interested in
/// @arg externalVertices: list of vertices which are external/rooted (vertex label AND whether we add L or L*)
template <typename T>
TestTemplatePureGauge<T>::TestTemplatePureGauge(const GraphContainer &container, const std::vector<ExternalPolyakovLoop> &externalVertices) :
    Container(container),
    ExternalVertices(externalVertices)
{
    if (this->ExternalVertices.size()>2)
        throw std::invalid_argument("ERROR: TestTemplatePureGauge requires externalVertices to be of size 2 or less!\n");

    if (std::find_if(this->ExternalVertices.begin(), this->ExternalVertices.end(), [this](const ExternalPolyakovLoop& x) { return (x.Label>this->Container.GetN() || x.Label<1); })!=this->ExternalVertices.end())
        throw std::invalid_argument("ERROR: TestTemplatePureGauge requires ExternalVertices to have labels between 1 and N!\n");

    for (int i=this->Container.GetNTimesNMinusOneDiv2()-1; i>=0; --i)
        if (this->Container.GetElementAdjacencyMatrix(this->Container.GetRowM(i),this->Container.GetColM(i)))
            this->Edges.push_back(UndirectedEdge(this->Container.GetRowM(i),this->Container.GetColM(i)));
}

/// binomial coefficient for doubles (recursive)
/// @arg n: first argument
/// @arg k: second argument
template <typename T>
double TestTemplatePureGauge<T>::BinomialCoefficient(int n, int k)
{
#ifdef DEBUG
    if (k > n)
        throw std::invalid_argument("BinomialCoefficient requires k <= n!\n");
    if (k < 0 || n <0)
        throw std::invalid_argument("BinomialCoefficient requires both n and k to be greater than zero!\n");
#endif
    if (k==0 || k==n)
        return 1;
    return (this->BinomialCoefficient(n-1, k-1) + this->BinomialCoefficient(n-1, k));
}

/// factorial for double (recursive)
/// @param n: size
template <typename T>
double TestTemplatePureGauge<T>::Factorial(int n)
{
#ifdef DEBUG
    if (n<0)
        throw std::invalid_argument("Factorial requires a non-negative argument!\n");
#endif
    return (n==0) || (n==1) ? 1 : n*this->Factorial(n-1);
}

/// generic template routine...
template <typename T>
T TestTemplatePureGauge<T>::SingleSiteIntegral(int n, int m)
{
    static_assert(sizeof(T)==0, "Only specializations of TestTemplatePureGauge::SingleSiteIntegral can be used!");
}

/// Single site integral: I_{n,m} = \int dU \chi(U)^n \chi(U^{\dagger})^m (see eq. (A.14) of J. Glesaaen's Thesis)
/// \chi(U) = e^{i\theta} + e^{i\phi} + e^{-i(\theta+\phi)}  = \chi(U^{\dagger})^*, for N_c=3
/// Glesaaen gives a recursion relation between I_{n,m} and J_{n,m}, where J_{n,m} is the integral over the angles (\theta, \phi) and does not include the Haar Measure
/// J. Scheunert has come up with a compact expression which does not use this recursion and this is what we implement below
/// I_{n,m} = \sum^{floor(n/3)}_{j=max{0, (n-m)/3} \frac {2 n! m!}{(n-j-(n-m)/3+1)!(n-j-(n-m)/3+2)!(2j-(n-m)/3)!} Binomial[2j-(n-m)/3,j] Binomial[3(n-j-(n-m)/3+1),n-3j]
/// of course if (n-m)%3 \neq 0 we get zero!
/// double specialization
template<>
double TestTemplatePureGauge<double>::SingleSiteIntegral(int n, int m)
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
        double bin1 = this->BinomialCoefficient(3*(n-j-nMinusMDiv3+1), n-3*j);
        double bin2 = this->BinomialCoefficient(2*j-nMinusMDiv3,j);
        double nFactorial = this->Factorial(n);
        double mFactorial = this->Factorial(m);
        double denom1 = this->Factorial(n-j-nMinusMDiv3+1);
        double denom2 = this->Factorial(n-j-nMinusMDiv3+2);
        double denom3 = this->Factorial(2*j-nMinusMDiv3);
        result += 2*(nFactorial/denom1)*(mFactorial/denom2)*(bin1*bin2/denom3);
    }
    return result;
}

/// GiNaC::numeric specialization (see description above)
template<>
GiNaC::numeric TestTemplatePureGauge<GiNaC::numeric>::SingleSiteIntegral(int n, int m)
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
        throw std::logic_error("TestTemplatePureGauge<GiNaC::numeric>::SingleSiteIntegral should give a rational number!\n");
    return result;
}

/// for a fixed set of directions for the edges of the graph, compute the weight \prod_i W_i, i=1,2,...,N_V where N_V is the number of vertices
/// @arg directedEdges: vector of length N_B which tells us "direction" of each link
template <typename T>
T TestTemplatePureGauge<T>::GetGraphWeightFixedBonds(const std::vector<bool>& directedEdges)
{
    auto siteCounts = this->GetSiteCountsForDirectedEdges(directedEdges);

    T result = 1;
    for (int i=0; i<siteCounts.size(); ++i)
    {
        auto siteContribution = this->SingleSiteIntegral(siteCounts[i].NbrIn, siteCounts[i].NbrOut);
#ifdef DEBUG
        std::cout << "Vertex " << i+1 << " has contribution " << siteContribution << "\n";
#endif
        result *= siteContribution;
    }
#ifdef DEBUG
    std::cout << "FINAL_RESULT: " << result << "\n";
#endif
    return result;
}

/// from a set of directed edges, get the counts of incoming and outcoming edges for each vertex
/// factor in external sites!
/// returns a vector of SiteCount objects of length N_V
/// @arg directedEdges: vector of bool variables of size N_B containing the edge direction
template <typename T>
inline std::vector<SiteCount> TestTemplatePureGauge<T>::GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges)
{
    if (directedEdges.size() != this->Edges.size())
        throw std::invalid_argument("PureGaugeWeight::GetSiteInformationGivenDirectedEdges requires directedEdges to be of the same size as Edges!\n");

    std::vector<SiteCount> result(this->Container.GetN(), {0,0});

    for (int i=0; i<directedEdges.size(); ++i)
    {
#ifdef DEBUG
        if (1 > this->Edges[i].FirstVertex || this->Edges[i].FirstVertex > this->Container.GetN() || 1 > this->Edges[i].SecondVertex || this->Edges[i].SecondVertex > this->Container.GetN())
            throw std::invalid_argument("PureGaugeWeight::GetSiteInformationGivenDirectedEdges requires Edges to connect vertices ranging from 1 to N!\n");
#endif
        if (directedEdges[i])
        {
            result[this->Edges[i].FirstVertex-1].NbrOut++;
            result[this->Edges[i].SecondVertex-1].NbrIn++;
#ifdef DEBUG
            std::cout << "EDGE " << i+1 << " from vertex " << this->Edges[i].FirstVertex << " to vertex " << this->Edges[i].SecondVertex << "\n";
#endif
        }
        else
        {
            result[this->Edges[i].FirstVertex-1].NbrIn++;
            result[this->Edges[i].SecondVertex-1].NbrOut++;
#ifdef DEBUG
            std::cout << "EDGE " << i+1 << " from vertex " << this->Edges[i].SecondVertex << " to vertex " << this->Edges[i].FirstVertex << "\n";
#endif
        }
    }

    /// loop over external vertices and add to in/out counts
    for (auto it=this->ExternalVertices.begin(); it!=this->ExternalVertices.end(); ++it)
    {
        if (it->Fundamental) /// add L
            result[it->Label-1].NbrIn++;
        else /// add L*
            result[it->Label-1].NbrOut++;
    }

#ifdef DEBUG
    for (int i=0; i<result.size(); ++i)
        std::cout << "Vertex " << i+1 << " has " << result[i].NbrIn << " vertices incoming and " << result[i].NbrOut << " outgoing!\n";
#endif
    return result;
}

/// recursive function which enumerates all possible combinations of edge directions
/// termination condition when nbrBondsRemaining = 0 and thus we compute contribution to TotalWeight and return
/// "true" means that bond is directed from FirstVertex to SecondVertex
/// "false" means that bond is directed to FirstVertex from SecondVertex
/// @arg tmp: vector containing a given combination of bond directions
/// @arg nbrBondsRemaining: N_B-k, where k is the size of tmp
template <typename T>
void TestTemplatePureGauge<T>::GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining)
{
    if (nbrBondsRemaining==0)
    {
        this->TotalWeight += this->GetGraphWeightFixedBonds(tmp);
        return;
    }

    /// true ("forward" i.e. first to second vertex)
    tmp.push_back(true);
    this->GetAllWeights(tmp, nbrBondsRemaining-1);
    tmp.pop_back();

    /// false ("backward" i.e. second to first vertex)
    tmp.push_back(false);
    this->GetAllWeights(tmp, nbrBondsRemaining-1);
    tmp.pop_back();

    return;
}

/// compute the weight of a given graph G
/// \phi(G) \equiv \int \prod_{x \in V(G)} dU(x) \prod_{(x,y) \in \sigma(E_G)} \left( L(x)L^*(y) + L^*(x)L(y) \right)
/// V(G) is the set of vertices, E_G is the set of edges, and \sigma is the edge to endpoint function.
/// This weight does not depend on the embedding of the graph inside the lattice
template <typename T>
inline T TestTemplatePureGauge<T>::Weight()
{
    std::vector<bool> tmp;
    this->TotalWeight = 0; /// set weight
    this->GetAllWeights(tmp, this->Edges.size());
    return this->TotalWeight;
}

// templated type can only be double or GiNaC::numeric
template class TestTemplatePureGauge<double>;
template class TestTemplatePureGauge<GiNaC::numeric>;
