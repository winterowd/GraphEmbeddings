#include "StaticQuarkWeight.h"

/// define static member variable
bool StaticQuarkWeight::StaticDeterminantPrecomputed;

/// constructor
/// @arg container: pointer to GraphContainer object describing the graph we are interested in
/// @arg externalVertices: list of vertices which are external/rooted (vertex label AND whether we add L or L*)
StaticQuarkWeight::StaticQuarkWeight(const GraphContainer &subgraphG, const GraphContainer &fullG, const std::vector<ExternalPolyakovLoop> &externalVertices) :
    Container(subgraphG),
    ExternalVertices(externalVertices)
{
    std::cout << "DEBUG_EXTERNALVERTICES: " << this->ExternalVertices.size() << "\n";
    if (this->ExternalVertices.size()>2)
        throw std::invalid_argument("ERROR: StaticQuarkWeight requires externalVertices to be of size 2 or less!\n");

    if (std::find_if(this->ExternalVertices.begin(), this->ExternalVertices.end(), [this](const ExternalPolyakovLoop& x) { return (x.Label>this->Container.GetN() || x.Label<1); })!=this->ExternalVertices.end())
        throw std::invalid_argument("ERROR: StaticQuarkWeight requires ExternalVertices to have labels between 1 and N!\n");

    this->DeltaNumberVertices = fullG.GetN()-this->Container.GetN();

    for (int i=this->Container.GetNTimesNMinusOneDiv2()-1; i>=0; --i)
        if (this->Container.GetElementAdjacencyMatrix(this->Container.GetRowM(i),this->Container.GetColM(i)))
            this->Edges.push_back(UndirectedEdge(this->Container.GetRowM(i),this->Container.GetColM(i)));
}

/// Single site integral: I_{n,m} = \int dU \chi(U)^n \chi(U^{\dagger})^m (see eq. (A.14) of J. Glesaaen's Thesis)
/// \chi(U) = e^{i\theta} + e^{i\phi} + e^{-i(\theta+\phi)}  = \chi(U^{\dagger})^*, for N_c=3
/// Glesaaen gives a recursion relation between I_{n,m} and J_{n,m}, where J_{n,m} is the integral over the angles (\theta, \phi) and does not include the Haar Measure
/// J. Scheunert has come up with a compact expression which does not use this recursion and this is what we implement below
/// I_{n,m} = \sum^{floor(n/3)}_{j=max{0, (n-m)/3} \frac {2 n! m!}{(n-j-(n-m)/3+1)!(n-j-(n-m)/3+2)!(2j-(n-m)/3)!} Binomial[2j-(n-m)/3,j] Binomial[3(n-j-(n-m)/3+1),n-3j]
/// of course if (n-m)%3 \neq 0 we get zero!
GiNaC::numeric StaticQuarkWeight::SingleSiteIntegral(int n, int m)
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
        throw std::logic_error("StaticQuarkWeight<GiNaC::numeric>::SingleSiteIntegral should give a rational number!\n");
    return result;
}

/// for a fixed set of directions for the edges of the graph, compute the weight \prod_i W_i, i=1,2,...,N_V where N_V is the number of vertices
/// @arg directedEdges: vector of length N_B which tells us "direction" of each link
GiNaC::ex StaticQuarkWeight::GetGraphWeightFixedBonds(const std::vector<bool>& directedEdges)
{
    auto siteCounts = this->GetSiteCountsForDirectedEdges(directedEdges);

    GiNaC::ex result = 1;
    for (int i=0; i<siteCounts.size(); ++i)
    {
        std::cout << "Vertex " << i+1 << " has " << siteCounts[i].NbrIn << " counts in and " << siteCounts[i].NbrOut << " out\n";
        auto siteContribution = this->ComputeSiteContribution(siteCounts[i].NbrIn, siteCounts[i].NbrOut);
#ifdef DEBUG
        std::cout << "Vertex " << i+1 << " has contribution " << siteContribution << "\n";
#endif
        std::cout << "Vertex " << i+1 << " has contribution " << siteContribution << "\n";
        result *= siteContribution;
    }
#ifdef DEBUG
    std::cout << "FINAL_RESULT: " << result << "\n";
#endif
    return result;
}

/// compute the contribution at a given site (vertex) including the static quark determinant with n incoming and m outgoing lines
/// @param nbrIn: number of "incoming" lines (powers of L)
/// @param nbrOut: number of "outoing" lines (powers of L^*)
GiNaC::ex StaticQuarkWeight::ComputeSiteContribution(int nbrIn, int nbrOut)
{
    GiNaC::ex result = GiNaC::numeric(0);
    for (int i=0; i<StaticQuarkWeight::StaticDeterminant.size(); ++i)
        result += StaticQuarkWeight::StaticDeterminant[i].first*this->SingleSiteIntegral(nbrIn+StaticQuarkWeight::StaticDeterminant[i].second.NbrIn, nbrOut+StaticQuarkWeight::StaticDeterminant[i].second.NbrOut);
    return result;
}

/// from a set of directed edges, get the counts of incoming and outcoming edges for each vertex
/// factor in external sites!
/// returns a vector of SiteCount objects of length N_V
/// @arg directedEdges: vector of bool variables of size N_B containing the edge direction
inline std::vector<SiteCount> StaticQuarkWeight::GetSiteCountsForDirectedEdges(const std::vector<bool> &directedEdges)
{
    if (directedEdges.size() != this->Edges.size())
        throw std::invalid_argument("StaticQuarkWeight::GetSiteInformationGivenDirectedEdges requires directedEdges to be of the same size as Edges!\n");

    std::vector<SiteCount> result(this->Container.GetN(), {0,0});

    for (int i=0; i<directedEdges.size(); ++i)
    {
#ifdef DEBUG
        if (1 > this->Edges[i].FirstVertex || this->Edges[i].FirstVertex > this->Container.GetN() || 1 > this->Edges[i].SecondVertex || this->Edges[i].SecondVertex > this->Container.GetN())
            throw std::invalid_argument("StaticQuarkWeight::GetSiteInformationGivenDirectedEdges requires Edges to connect vertices ranging from 1 to N!\n");
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
void StaticQuarkWeight::GetAllWeights(std::vector<bool>& tmp, int nbrBondsRemaining)
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

/// compute the weight of a given graph G: \phi(g) = z^{|V(G)|-|V(g)|}_0 \tilde{\phi}(g)
/// \tilde{\phi}(g) \equiv \int \prod_{x \in V(g)} dU(x) \prod_{(x,y) \in \sigma(E_g)} \left( L(x)L^*(y) + L^*(x)L(y) \right)
/// z_0 \equiv \int dL_x det_x Q_stat, (integrated static determinant)
/// V(G) is the set of vertices, E_G is the set of edges, and \sigma is the edge to endpoint function.
/// This weight does not depend on the embedding of the graph inside the lattice
GiNaC::ex StaticQuarkWeight::Weight()
{
     /// check if integrated and unintegrated static determinant has been computed
    if (!StaticQuarkWeight::StaticDeterminantPrecomputed)
    {
        StaticQuarkWeight::StaticDeterminant = ComputeStaticDeterminant();
        SingleSiteIntegratedStaticDeterminant = ComputeSingleSiteIntegratedStaticDeterminant();
        StaticQuarkWeight::StaticDeterminantPrecomputed = true;
    }

    std::vector<bool> tmp;
    this->TotalWeight = GiNaC::numeric(0); /// initialize weight
    this->GetAllWeights(tmp, this->Edges.size()); /// \tilde{\phi}(g)
    this->TotalWeight *= GiNaC::pow(StaticQuarkWeight::SingleSiteIntegratedStaticDeterminant, this->DeltaNumberVertices); /// z^{|V(G)|-|V(g)|}_0
    return this->TotalWeight;
}

/// compute the single site static determinant (hard-coded)
/// \det_x Q_stat = (1 + h_1 L + h^2_1 L^* + h^3_1)^2 (1 + \bar{h}_1 L + \bar{h}^2_1 L^* + \bar{h}^3_1)^2
std::vector<std::pair<GiNaC::ex, SiteCount>> ComputeStaticDeterminant()
{
    /// for debugging purposes
    /*std::vector<std::pair<double,double>> testPoints;
    testPoints.push_back(std::pair<double,double>(0.12398,8.98));
    testPoints.push_back(std::pair<double,double>(-0.12398,8.98));
    testPoints.push_back(std::pair<double,double>(0.12398,-8.98));
    testPoints.push_back(std::pair<double,double>(0.549,0.000023));
    testPoints.push_back(std::pair<double,double>(0.00823,0.010023));*/

    std::vector<std::pair<GiNaC::ex, SiteCount>> result;
    GiNaC::symbol h1 = AuxiliaryRoutinesForGinac::GetSymbol(0, false);
    GiNaC::symbol hBar1 = AuxiliaryRoutinesForGinac::GetSymbol(1, false);
    std::vector<GiNaC::ex> powersH1, powersHBar1;
    for (int i=0; i<6; ++i)
    {
        powersH1.push_back(GiNaC::pow(h1,i+1));
        powersHBar1.push_back(GiNaC::pow(hBar1,i+1));
    }
    GiNaC::ex temp = GiNaC::numeric(1)+GiNaC::numeric(2)*powersH1[2]+powersH1[5]+powersHBar1[5]+GiNaC::numeric(2)*powersHBar1[2]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[2]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[2]+GiNaC::numeric(2)*powersH1[2]*powersHBar1[5]+powersH1[5]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (0,0)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{0,0}));
    temp = GiNaC::numeric(2)*powersH1[0]+GiNaC::numeric(2)*powersH1[3]+GiNaC::numeric(2)*powersHBar1[1]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[2]+GiNaC::numeric(2)*powersHBar1[4]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[4]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[4]+GiNaC::numeric(2)*powersH1[0]*powersHBar1[5]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (1,0)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{1,0}));
    temp = powersH1[1]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[1]*powersHBar1[2]+powersHBar1[3]+GiNaC::numeric(2)*powersH1[2]*powersHBar1[3]+powersH1[5]*powersHBar1[3]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[4]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[4]+powersH1[1]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (2,0)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{2,0}));
    temp = GiNaC::numeric(2)*powersH1[1]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[0]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[1]*powersHBar1[4];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (3,0)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{3,0}));
    result.push_back(std::pair<GiNaC::ex, SiteCount>(powersH1[1]*powersHBar1[3], SiteCount{4,0}));
    temp = GiNaC::numeric(2)*powersHBar1[0]+GiNaC::numeric(2)*powersH1[1]+GiNaC::numeric(2)*powersH1[4]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[0]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[0]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[2]+GiNaC::numeric(2)*powersHBar1[3]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[1]*powersHBar1[5]+GiNaC::numeric(2)*powersH1[4]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (0,1)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{0,1}));
    temp = powersHBar1[1]+powersH1[3]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[0]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[0]+GiNaC::numeric(2)*powersH1[2]*powersHBar1[1]+powersH1[5]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[3]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[3]+powersH1[3]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (0,2)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{0,2}));
    temp = GiNaC::numeric(2)*powersH1[3]*powersHBar1[0]+GiNaC::numeric(2)*powersH1[1]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[4]*powersHBar1[1];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (0,3)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{0,3}));
    result.push_back(std::pair<GiNaC::ex, SiteCount>(powersH1[3]*powersHBar1[1], SiteCount{0,4}));
    temp = GiNaC::numeric(2)*powersH1[2]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[0]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[0]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[1]+GiNaC::numeric(2)*powersHBar1[2]+GiNaC::numeric(8)*powersH1[2]*powersHBar1[2]+GiNaC::numeric(2)*powersH1[5]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[3]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[3]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[4]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[4]+GiNaC::numeric(2)*powersH1[2]*powersHBar1[5];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (1,1)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{1,1}));
    temp = powersH1[1]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[2]+powersH1[3]*powersHBar1[3];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (2,2)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{2,2}));
    temp = GiNaC::numeric(4)*powersH1[2]*powersHBar1[0]+GiNaC::numeric(2)*powersH1[0]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[4];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (1,2)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{1,2}));
    temp = GiNaC::numeric(2)*powersH1[1]*powersHBar1[0]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[1]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[3]*powersHBar1[2]+GiNaC::numeric(4)*powersH1[1]*powersHBar1[3]+GiNaC::numeric(2)*powersH1[4]*powersHBar1[3]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[4];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (2,1)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{2,1}));
    temp = GiNaC::numeric(2)*powersH1[2]*powersHBar1[1]+GiNaC::numeric(2)*powersH1[3]*powersHBar1[2];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (1,3)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{1,3}));
    temp = GiNaC::numeric(2)*powersH1[1]*powersHBar1[2]+GiNaC::numeric(2)*powersH1[2]*powersHBar1[3];
    /*std::cout << "DEBUG_STATIC_DETERMINANT: (3,1)\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << temp.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    result.push_back(std::pair<GiNaC::ex, SiteCount>(temp, SiteCount{3,1}));
    /*for (int i=0; i<result.size(); ++i)
        std::cout << result[i].first << " (" << result[i].second.NbrIn << "," << result[i].second.NbrOut << ")\n";*/
    return result;
}

/// compute the single site integrated static determinant (hard-coded)
/// \int dL_x det_x Q_stat
GiNaC::ex ComputeSingleSiteIntegratedStaticDeterminant()
{
    /// for debugging purposes
    /*std::vector<std::pair<double,double>> testPoints;
    testPoints.push_back(std::pair<double,double>(0.12398,8.98));
    testPoints.push_back(std::pair<double,double>(-0.12398,8.98));
    testPoints.push_back(std::pair<double,double>(0.12398,-8.98));
    testPoints.push_back(std::pair<double,double>(0.549,0.000023));
    testPoints.push_back(std::pair<double,double>(0.00823,0.010023));*/

    GiNaC::symbol h1 = AuxiliaryRoutinesForGinac::GetSymbol(0, false);
    GiNaC::symbol hBar1 = AuxiliaryRoutinesForGinac::GetSymbol(1, false);
    std::vector<GiNaC::ex> powersH1, powersHBar1;
    for (int i=0; i<6; ++i)
    {
        powersH1.push_back(GiNaC::pow(h1,i+1));
        powersHBar1.push_back(GiNaC::pow(hBar1,i+1));
    }
    GiNaC::ex result = 1+GiNaC::numeric(4)*powersH1[2]+powersH1[5]+GiNaC::numeric(4)*powersHBar1[2]+powersHBar1[5]+GiNaC::numeric(4)*powersH1[0]*powersHBar1[0]+GiNaC::numeric(6)*powersH1[3]*powersHBar1[0]+GiNaC::numeric(10)*powersH1[1]*powersHBar1[1]+GiNaC::numeric(6)*powersH1[4]*powersHBar1[1]+GiNaC::numeric(20)*powersH1[2]*powersHBar1[2];
    result += GiNaC::numeric(4)*powersH1[5]*powersHBar1[2]+GiNaC::numeric(6)*powersH1[0]*powersHBar1[3]+GiNaC::numeric(10)*powersH1[3]*powersHBar1[3]+GiNaC::numeric(6)*powersH1[1]*powersHBar1[4]+GiNaC::numeric(4)*powersH1[4]*powersHBar1[4]+GiNaC::numeric(4)*powersH1[2]*powersHBar1[5]+powersH1[5]*powersHBar1[5];
    /*std::cout << "DEBUG_INTEGRATED_STATIC_DETERMINANT:\n";
    for (int i=0; i<testPoints.size(); ++i)
        std::cout << result.subs(h1 == testPoints[i].first).subs(hBar1 == testPoints[i].second) << " ";
    std::cout << "\n";*/
    return result;
}

/// static member variables
std::vector<std::pair<GiNaC::ex, SiteCount>> StaticQuarkWeight::StaticDeterminant;
GiNaC::ex StaticQuarkWeight::SingleSiteIntegratedStaticDeterminant;
