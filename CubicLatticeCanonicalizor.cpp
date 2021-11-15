#include "CubicLatticeCanonicalizor.h"

/// constructor for CubicLatticeCanonicalizor
/// @input container: pointer to GraphContainer object which describes connectivity etc
/// @input lattice: pointer to CubicLattice object (not AbstractLattice)
/// @input embedList: copy of VertexEmbedList so we can know where the vertices are on the lattice
CubicLatticeCanonicalizor::CubicLatticeCanonicalizor(GraphContainer *container, CubicLattice *lattice, const VertexEmbedList &embedList) :
    Container(container),
    Lattice(lattice),
    VerticesCartesian(container->GetN(), {0,0,0}),
    VerticesCartesianCOM(container->GetN(), {0,0,0}),
    COM({0,0,0}),
    OriginalList(embedList)
{
    if (embedList.GetSize()!=container->GetN())
        throw std::invalid_argument("CubicLatticeCanonicalizor requires embedList to be of the same size as container!\n");

    std::vector<unsigned int> indices(3, 0);
    for (auto it=embedList.begin(); it!=embedList.end(); ++it)
    {
        this->Lattice->GetSiteCoordinates(it->Index, indices);
        /// NOTE: assume that vertex numbers start at 1 in VertexEmbedList so need to shift!
        this->VerticesCartesian[it->Number-1] = CartesianCoords{double(indices[0]), double(indices[1]), double(indices[2])};
        /// TODO: check if graph is rotoed from VertexEmbedList and then add appropriate color to CartesianCoords (prepend)
    }

}

/// a += b
void CubicLatticeCanonicalizor::AccumulateCartesianVector(CartesianCoords& a, const CartesianCoords& b)
{
    for (int i=0; i<3; ++i)
        a[i] += b[i];
}

/// a -= b
void CubicLatticeCanonicalizor::AccumulateCartesianVectorNeg(CartesianCoords& a, const CartesianCoords& b)
{
    for (int i=0; i<3; ++i)
        a[i] -= b[i];
}

/// a *= s
void CubicLatticeCanonicalizor::ScalarMultCartesianVector(CartesianCoords& a, double s)
{
    for (int i=0; i<3; ++i)
        a[i] *= s;
}

/// norm squared of a vector
double CubicLatticeCanonicalizor::NormSqCartesianVector(const CartesianCoords& a)
{
    double result=0;
    for (int i=0; i<3; ++i)
        result += a[i]*a[i];
    return result;
}

/// norm a vector
double CubicLatticeCanonicalizor::NormCartesianVector(const CartesianCoords& a)
{
    return std::sqrt(this->NormSqCartesianVector(a));
}

/// compute the COM vector and then shift vertices such that COM lies at the origin
/// COM vector = (\sum_v_i \vec{v}_i + \frac{1}{2} \sum_E_i ( \vec{i}(E_i) + \vec{f}(E_i) ) )/(|V| + |E|)
/// result put into VerticesCartesianCOM
void CubicLatticeCanonicalizor::ComputeCOMandShift()
{

    this->VerticesCartesianCOM = this->VerticesCartesian; /// copy vertices into COM

    this->COM = {0,0,0}; /// clear COM vector

    /// compute the COM vector
    CartesianCoords temp{0,0,0};
    for (int i=this->Container->GetNTimesNMinusOneDiv2()-1; i>=0; --i)
    {
        if (this->Container->GetElementAdjacencyMatrix(this->Container->GetRowM(i),this->Container->GetColM(i)))
        {
            /// NOTE: vertex numbers start at 1 in GraphContainer so need to shift index!
            this->AccumulateCartesianVector(temp,this->VerticesCartesian[this->Container->GetRowM(i)-1]);
            this->AccumulateCartesianVector(temp,this->VerticesCartesian[this->Container->GetColM(i)-1]);
        }
    }
    this->ScalarMultCartesianVector(temp, 0.5);
    this->AccumulateCartesianVector(this->COM, temp);
    for (int i=0; i<this->VerticesCartesian.size(); ++i)
        this->AccumulateCartesianVector(this->COM, this->VerticesCartesian[i]);
    this->ScalarMultCartesianVector(this->COM, 1./(this->VerticesCartesian.size()+this->Container->GetL()));

    /// subtract the vector corresponding to the the COM in the old coordinate system
    for (int i=0; i<this->VerticesCartesian.size(); ++i)
        this->AccumulateCartesianVectorNeg(this->VerticesCartesianCOM[i], this->COM);
}

/// rotation of vector x about axis (unit-length) passing through point p
CubicLatticeCanonicalizor::CartesianCoords CubicLatticeCanonicalizor::RotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords& p, const CartesianCoords &x)
{
    if (std::fabs(this->NormCartesianVector(axis)-1.)>std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("RotationTransform requires axis to be of unit length!\n");

    CartesianCoords temp(x);
    CartesianCoords result{0,0,0};

    /// translate by -p
    this->AccumulateCartesianVectorNeg(temp, p);

    /// rotate by theta about axis
    double cosTheta = std::cos(theta);
    double oneMinusCosTheta = 1-cosTheta;
    double sinTheta = std::sin(theta);

    /// hard-code R times the vector
    result[0] = (cosTheta+axis[0]*axis[0]*oneMinusCosTheta)*temp[0]+(axis[0]*axis[1]*oneMinusCosTheta-axis[2]*sinTheta)*temp[1]+(axis[0]*axis[2]*oneMinusCosTheta+axis[1]*sinTheta)*temp[2];
    result[1] = (axis[0]*axis[1]*oneMinusCosTheta+axis[2]*sinTheta)*temp[0]+(cosTheta+axis[1]*axis[1]*oneMinusCosTheta)*temp[1]+(axis[1]*axis[2]*oneMinusCosTheta-axis[0]*sinTheta)*temp[2];
    result[2] = (axis[0]*axis[2]*oneMinusCosTheta-axis[1]*sinTheta)*temp[0]+(axis[2]*axis[1]*oneMinusCosTheta+axis[0]*sinTheta)*temp[1]+(cosTheta+axis[2]*axis[2]*oneMinusCosTheta)*temp[2];

    /// translate by p
    this->AccumulateCartesianVector(temp, p);

    return result;
}

/// rotation of vector x about axis (unit-length) passing through the origin
CubicLatticeCanonicalizor::CartesianCoords CubicLatticeCanonicalizor::RotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords &x)
{
    if (std::fabs(this->NormCartesianVector(axis)-1.)>std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("RotationTransform requires axis to be of unit length!\n");

    CartesianCoords result{0,0,0};

    /// rotate by theta about axis
    double cosTheta = std::cos(theta);
    double oneMinusCosTheta = 1-cosTheta;
    double sinTheta = std::sin(theta);

    /// hard-code R times the vector
    result[0] = (cosTheta+axis[0]*axis[0]*oneMinusCosTheta)*x[0]+(axis[0]*axis[1]*oneMinusCosTheta-axis[2]*sinTheta)*x[1]+(axis[0]*axis[2]*oneMinusCosTheta+axis[1]*sinTheta)*x[2];
    result[1] = (axis[0]*axis[1]*oneMinusCosTheta+axis[2]*sinTheta)*x[0]+(cosTheta+axis[1]*axis[1]*oneMinusCosTheta)*x[1]+(axis[1]*axis[2]*oneMinusCosTheta-axis[0]*sinTheta)*x[2];
    result[2] = (axis[0]*axis[2]*oneMinusCosTheta-axis[1]*sinTheta)*x[0]+(axis[2]*axis[1]*oneMinusCosTheta+axis[0]*sinTheta)*x[1]+(cosTheta+axis[2]*axis[2]*oneMinusCosTheta)*x[2];

    return result;
}

/// inversion about the point p
/// inv(x) = 2p - x
CubicLatticeCanonicalizor::CartesianCoords CubicLatticeCanonicalizor::Inversion(const CartesianCoords& p, const CartesianCoords& x)
{
    CartesianCoords result(p);
    this->ScalarMultCartesianVector(result,2.);
    this->AccumulateCartesianVectorNeg(result,x);
    return result;
}

/// rotation followed by inversion
/// inv(R(x))
CubicLatticeCanonicalizor::CartesianCoords CubicLatticeCanonicalizor::InversionComposedWithRotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords& p, const CartesianCoords& x)
{
    return this->Inversion(p,this->RotationTransform(theta, axis, p, x));
}

/// apply rotations to a Cartesian vector about four-fold axes which pass through the centers of opposing faces
std::vector<std::vector<CubicLatticeCanonicalizor::CartesianCoords>> CubicLatticeCanonicalizor::ApplyRotationsFourFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices)
{
    std::vector<std::vector<CartesianCoords>> result;

    std::vector<CartesianCoords> axes{{1,0,0}, {0,1,0}, {0,0,1}}; /// axes
    std::vector<double> angles{-M_PI_2,M_PI_2,M_PI}; /// angles (4-fold)

    for (auto axes_it=axes.begin(); axes_it!=axes.end(); ++axes_it) ///  loop over axes
    {
        for (auto angles_it=angles.begin(); angles_it!=angles.end(); ++angles_it) /// loop over angles
        {
            std::vector<CartesianCoords> rotate;
            for (int i=0; i<vertices.size(); ++i)
                rotate.push_back(this->RotationTransform(*angles_it, *axes_it, p, vertices[i])); /// rotation
            result.push_back(rotate);
            std::vector<CartesianCoords> rotatePlusInversion;
            for (int i=0; i<vertices.size(); ++i)
                rotatePlusInversion.push_back(this->InversionComposedWithRotationTransform(*angles_it, *axes_it, p, vertices[i])); /// rotation followed by inversion
            result.push_back(rotatePlusInversion);
        }
    }

    return result;
}

/// apply rotations to a Cartesian vector about two-fold axes which pass bewteen the midpoints of diagonally opposing edges
std::vector<std::vector<CubicLatticeCanonicalizor::CartesianCoords>> CubicLatticeCanonicalizor::ApplyRotationsTwoFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices)
{
    std::vector<std::vector<CartesianCoords>> result;

    double inverseSqrtTwo = 1./std::sqrt(2);
    std::vector<CartesianCoords> axes{{0,inverseSqrtTwo,inverseSqrtTwo}, {0,-inverseSqrtTwo,inverseSqrtTwo}, {inverseSqrtTwo,0,inverseSqrtTwo},
                                      {inverseSqrtTwo,0,-inverseSqrtTwo}, {inverseSqrtTwo,inverseSqrtTwo,0}, {-inverseSqrtTwo,inverseSqrtTwo,0}}; /// axes
    double theta = M_PI; /// (2-fold)
    for (auto axes_it=axes.begin(); axes_it!=axes.end(); ++axes_it) ///  loop over axes
    {
        std::vector<CartesianCoords> rotate;
        for (int i=0; i<vertices.size(); ++i)
            rotate.push_back(this->RotationTransform(theta, *axes_it, p, vertices[i])); /// rotation
        result.push_back(rotate);
        std::vector<CartesianCoords> rotatePlusInversion;
        for (int i=0; i<vertices.size(); ++i)
            rotatePlusInversion.push_back(this->InversionComposedWithRotationTransform(theta, *axes_it, p, vertices[i])); /// rotation followed by inversion
        result.push_back(rotatePlusInversion);
    }
    return result;
}

/// apply rotations to a Cartesian vector about three-fold axes which pass betwen diagonally opposing vertices
std::vector<std::vector<CubicLatticeCanonicalizor::CartesianCoords>> CubicLatticeCanonicalizor::ApplyRotationsThreeFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices)
{
    std::vector<std::vector<CartesianCoords>> result;

    double inverseSqrtThree = 1./std::sqrt(3);
    std::vector<CartesianCoords> axes{{inverseSqrtThree,inverseSqrtThree,inverseSqrtThree}, {inverseSqrtThree,-inverseSqrtThree,inverseSqrtThree},
                                      {-inverseSqrtThree,inverseSqrtThree,inverseSqrtThree}, {-inverseSqrtThree,-inverseSqrtThree,inverseSqrtThree}}; /// axes
    std::vector<double> angles{-2*M_PI/3.,2.*M_PI/3.}; /// angles (3-fold)

    for (auto axes_it=axes.begin(); axes_it!=axes.end(); ++axes_it) ///  loop over axes
    {
        for (auto angles_it=angles.begin(); angles_it!=angles.end(); ++angles_it) /// loop over angles
        {
            std::vector<CartesianCoords> rotate;
            for (int i=0; i<vertices.size(); ++i)
                rotate.push_back(this->RotationTransform(*angles_it, *axes_it, p, vertices[i])); /// rotation
            result.push_back(rotate);
            std::vector<CartesianCoords> rotatePlusInversion;
            for (int i=0; i<vertices.size(); ++i)
                rotatePlusInversion.push_back(this->InversionComposedWithRotationTransform(*angles_it, *axes_it, p, vertices[i])); /// rotation followed by inversion
            result.push_back(rotatePlusInversion);
        }
    }
    return result;
}

void CubicLatticeCanonicalizor::GenerateAllPermutationsWithoutRepeats(std::vector<std::vector<double>> &lists, const std::vector<double>& s, std::vector<double>& temp)
{
    if (temp.size() == 3) /// we want only three coordinates
    {
        /// check if duplicates
        if (std::fabs(temp[0]) != std::fabs(temp[1]) && std::fabs(temp[0]) != std::fabs(temp[2]) && std::fabs(temp[1]) != std::fabs(temp[2]))
            lists.push_back(temp);
        return;
    }
    for (int i=0; i<s.size(); ++i)
    {
        temp.push_back(s[i]);
        GenerateAllPermutationsWithoutRepeats(lists, s, temp);
        temp.pop_back();
    }
}

/// generate all permutations which correspond to the 48 operations of the octahedral group
/// given a vector r=(x,y,z), we take all permutations of (\pm, \pm y, \pm z), where the signs can be taken independently of one another
/// @input lists: list of all the permutations for the graph (return)
/// @input s: the graph vertices
void CubicLatticeCanonicalizor::GenerateAllPermutationsWithoutRepeats(std::vector<CartesianCoords> &lists, const std::vector<double>& s)
{
    if (s.size()!=6)
        throw std::invalid_argument("GenerateAllPermutationsWithoutRepeats expects s to be of size 6!\n");

    lists.push_back(CartesianCoords{s[0],s[1],s[2]}); lists.push_back(CartesianCoords{s[0],s[2],s[1]}); lists.push_back(CartesianCoords{s[1],s[0],s[2]}); lists.push_back(CartesianCoords{s[1],s[2],s[0]}); lists.push_back(CartesianCoords{s[2],s[0],s[1]}); lists.push_back(CartesianCoords{s[2],s[1],s[0]});
    lists.push_back(CartesianCoords{s[3],s[1],s[2]}); lists.push_back(CartesianCoords{s[3],s[2],s[1]}); lists.push_back(CartesianCoords{s[1],s[3],s[2]}); lists.push_back(CartesianCoords{s[1],s[2],s[3]}); lists.push_back(CartesianCoords{s[2],s[3],s[1]}); lists.push_back(CartesianCoords{s[2],s[1],s[3]});
    lists.push_back(CartesianCoords{s[0],s[4],s[2]}); lists.push_back(CartesianCoords{s[0],s[2],s[4]}); lists.push_back(CartesianCoords{s[4],s[0],s[2]}); lists.push_back(CartesianCoords{s[4],s[2],s[0]}); lists.push_back(CartesianCoords{s[2],s[0],s[4]}); lists.push_back(CartesianCoords{s[2],s[4],s[0]});
    lists.push_back(CartesianCoords{s[0],s[1],s[5]}); lists.push_back(CartesianCoords{s[0],s[5],s[1]}); lists.push_back(CartesianCoords{s[1],s[0],s[5]}); lists.push_back(CartesianCoords{s[1],s[5],s[0]}); lists.push_back(CartesianCoords{s[5],s[0],s[1]}); lists.push_back(CartesianCoords{s[5],s[1],s[0]});
    lists.push_back(CartesianCoords{s[3],s[4],s[2]}); lists.push_back(CartesianCoords{s[3],s[2],s[4]}); lists.push_back(CartesianCoords{s[4],s[3],s[2]}); lists.push_back(CartesianCoords{s[4],s[2],s[3]}); lists.push_back(CartesianCoords{s[2],s[3],s[4]}); lists.push_back(CartesianCoords{s[2],s[4],s[3]});
    lists.push_back(CartesianCoords{s[3],s[1],s[5]}); lists.push_back(CartesianCoords{s[3],s[5],s[1]}); lists.push_back(CartesianCoords{s[1],s[3],s[5]}); lists.push_back(CartesianCoords{s[1],s[5],s[3]}); lists.push_back(CartesianCoords{s[5],s[3],s[1]}); lists.push_back(CartesianCoords{s[5],s[1],s[3]});
    lists.push_back(CartesianCoords{s[0],s[4],s[5]}); lists.push_back(CartesianCoords{s[0],s[5],s[4]}); lists.push_back(CartesianCoords{s[4],s[0],s[5]}); lists.push_back(CartesianCoords{s[4],s[5],s[0]}); lists.push_back(CartesianCoords{s[5],s[0],s[4]}); lists.push_back(CartesianCoords{s[5],s[4],s[0]});
    lists.push_back(CartesianCoords{s[3],s[4],s[5]}); lists.push_back(CartesianCoords{s[3],s[5],s[4]}); lists.push_back(CartesianCoords{s[4],s[3],s[5]}); lists.push_back(CartesianCoords{s[4],s[5],s[3]}); lists.push_back(CartesianCoords{s[5],s[3],s[4]}); lists.push_back(CartesianCoords{s[5],s[4],s[3]});
}

/// generate the appropriate permutations on all vertices of the graph
/// return a vector of vectors of CartesianCoords objects where the first index corresponds to the permutation number and the second corresponds to the vertex of the graph
/// @input vertices: vertices of the graph
std::vector<std::vector<CubicLatticeCanonicalizor::CartesianCoords>> CubicLatticeCanonicalizor::PermutationsOnAllVertices(const std::vector<CartesianCoords>& vertices)
{
    std::vector<std::vector<CartesianCoords>> result(48 /* number of permutations (hard-coded) */, std::vector<CartesianCoords>(vertices.size()));
    for (int i=0; i<vertices.size(); ++i) /// loop over vertices
    {
        std::vector<double> coordsExtended(6,-1);
        std::copy(vertices[i].begin(), vertices[i].end(), coordsExtended.begin());
        for (int j=0; j<vertices[i].size(); ++j)
            coordsExtended[j+3] = -vertices[i][j];
        std::vector<CartesianCoords> temp;
        this->GenerateAllPermutationsWithoutRepeats(temp, coordsExtended); /// permutations of coordinates and their negative
        for (int j=0; j<temp.size(); ++j)
            result[j][i] = temp[j];
    }

    return result;
}

/// old version of the canonicalization where we use the explicit rotations about the symmetry axes of the cube
/// lexicographical order determines the canonical graph
VertexEmbedList CubicLatticeCanonicalizor::GetCanonicalOld()
{
    CartesianCoords origin{0,0,0};

    this->ComputeCOMandShift();

    auto resultTwoFold = this->ApplyRotationsTwoFoldAxisOnGraph(origin, this->VerticesCartesianCOM);
    std::vector<CartesianCoords> canonical = *std::min_element(resultTwoFold.begin(), resultTwoFold.end());

    auto resultThreeFold = this->ApplyRotationsThreeFoldAxisOnGraph(origin, this->VerticesCartesianCOM);
    if (*std::min_element(resultThreeFold.begin(), resultThreeFold.end()) < canonical)
        canonical = *std::min_element(resultThreeFold.begin(), resultThreeFold.end());

    auto resultFourFold = this->ApplyRotationsFourFoldAxisOnGraph(origin, this->VerticesCartesianCOM);
    if (*std::min_element(resultFourFold.begin(), resultFourFold.end()) < canonical)
        canonical = *std::min_element(resultFourFold.begin(), resultFourFold.end());

    /// we now have the canonical graph. now need to get position of smallest lexicographical vertex
    CartesianCoords minVector = *std::min_element(canonical.begin(), canonical.end());

    CartesianCoords midpoint{this->Lattice->GetN()/2., this->Lattice->GetN()/2., this->Lattice->GetN()/2.};

    /// translate all members of canonical graph by -minVector (smallest lexicographical vertex now at the origin)
    for (int i=0; i<canonical.size(); ++i)
    {
        this->AccumulateCartesianVectorNeg(canonical[i],minVector);
        this->AccumulateCartesianVector(canonical[i], midpoint);
    }

    return this->ConvertCartesianVectorToVertexEmbedList(canonical);
}

/// new method where we use permutations to produce the canonical graph
/// lexicographical order also used here
VertexEmbedList CubicLatticeCanonicalizor::GetCanonical()
{
    this->ComputeCOMandShift(); /// shift to COM

    auto resultPermuations = this->PermutationsOnAllVertices(this->VerticesCartesianCOM); /// apply permutations (should be equivalent to rotations)

    std::vector<CartesianCoords> canonical = *std::min_element(resultPermuations.begin(), resultPermuations.end()); /// find canonical embedding

    /// we now have the canonical graph. now need to get position of smallest lexicographical vertex
    CartesianCoords minVector = *std::min_element(canonical.begin(), canonical.end());

    CartesianCoords midpoint{this->Lattice->GetN()/2., this->Lattice->GetN()/2., this->Lattice->GetN()/2.}; /// midpoint of lattice

    /// translate all members of canonical graph by -minVector (smallest lexicographical vertex now at the origin) and midpoint of lattice (N/2,N/2,N/2)
    for (int i=0; i<canonical.size(); ++i)
    {
        this->AccumulateCartesianVectorNeg(canonical[i],minVector);
        this->AccumulateCartesianVector(canonical[i], midpoint);
    }

    return this->ConvertCartesianVectorToVertexEmbedList(canonical);
}

/// convert the vector of CartesianCoordinates into a VertexEmbedList object
VertexEmbedList CubicLatticeCanonicalizor::ConvertCartesianVectorToVertexEmbedList(const std::vector<CartesianCoords>& vec)
{
    VertexEmbedList result = this->OriginalList.IsTwoPointFunction() ? VertexEmbedList(this->OriginalList.GetMaxLength(), this->OriginalList.GetCorrelatorDistance()) : VertexEmbedList(this->OriginalList.GetMaxLength());

    if (this->OriginalList.IsTwoPointFunction()) /// correlator
    {
        std::vector<int> fixedVertexNumbers{this->OriginalList.GetFixedVertex(0).Number,this->OriginalList.GetFixedVertex(1).Number};
        for (int i=0; i<vec.size(); ++i)
        {
            std::vector<unsigned int> indicesFromDouble(3);
            for (int j=0; j<3; ++j)
            {
                if (vec[i][j]<0)
                    throw std::invalid_argument("CubicLatticeCanonicalizor::ConvertCartesianVectorToVertexEmbedList converting from double to unsigned int using imaginary values!\n");
                indicesFromDouble[j] = vec[i][j];
            }
            auto index = this->Lattice->GetSiteIndex(indicesFromDouble);
            if (i==fixedVertexNumbers[0])
                result.AddFixedVertexEmbed(0,i,index);
            else if (i==fixedVertexNumbers[1])
                result.AddFixedVertexEmbed(1,i,index);
            else
                result.AddVertexEmbed(i,index);
        }
    }
    else /// otherwise
    {
        for (int i=0; i<vec.size(); ++i)
        {
            std::vector<unsigned int> indicesFromDouble(3);
            for (int j=0; j<3; ++j)
            {
                if (vec[i][j]<0)
                    throw std::invalid_argument("CubicLatticeCanonicalizor::GetCanonical converting from double to unsigned int using imaginary values!\n");
                indicesFromDouble[j] = vec[i][j];
            }
            auto index = this->Lattice->GetSiteIndex(indicesFromDouble);
            result.AddVertexEmbed(i,index);
        }
    }
    return result;
}

void CubicLatticeCanonicalizor::PrintVerticesCartesian()
{
    std::cout << "**** VerticesCartesian ****\n";
    for (int i=0; i<this->VerticesCartesian.size(); ++i)
        std::cout << this->VerticesCartesian[i] << "\n";
}

void CubicLatticeCanonicalizor::PrintVerticesCartesianCOM()
{
    std::cout << "**** VerticesCartesianCOM ****\n";
    for (int i=0; i<this->VerticesCartesianCOM.size(); ++i)
        std::cout << this->VerticesCartesianCOM[i] << "\n";
}

/// lexicographical order
bool operator<(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs)
{
    for (int i=0; i<3; ++i)
        if (std::fabs(lhs[i]-rhs[i])>std::numeric_limits<double>::epsilon())
            return lhs[i]<rhs[i];
    return false;
}

/// equality
bool operator==(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs)
{
    for (int i=0; i<3; ++i)
        if (std::fabs(lhs[i]-rhs[i])>std::numeric_limits<double>::epsilon())
            return false;
    return true;
}

bool operator!=(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs)
{
    return !(lhs==rhs);
}

bool operator<(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs)
{
    if (lhs.size()!=rhs.size())
        return lhs.size()<rhs.size();
    for (int i=0; i<lhs.size(); ++i)
        if (lhs[i]!=rhs[i])
            return lhs[i]<rhs[i];
    return false;
}

bool operator==(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs)
{
    if (lhs.size()!=rhs.size())
        return false;
    for (int i=0; i<lhs.size(); ++i)
        if (lhs[i]!=rhs[i])
            return false;
    return true;
}

bool operator!=(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs)
{
    return !(lhs==rhs);
}

/// use lexicographical order for a pair which contains the Coordinates and the vertex number
/// need to keep the vertex number of each embedded vertex but VertexEmbed object sorts in lexicographical order with vertex number taking the first entry and site index the second
bool operator<(const std::pair<CubicLatticeCanonicalizor::CartesianCoords,int>& lhs, const std::pair<CubicLatticeCanonicalizor::CartesianCoords,int>& rhs)
{
    return (lhs.first<rhs.first);
}

std::ostream& operator<<(std::ostream& os, const CubicLatticeCanonicalizor::CartesianCoords &list)
{
    std::cout << "(" << list[0] << "," << list[1] << "," << list[2] << ")";
    return os;
}
