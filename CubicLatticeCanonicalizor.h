#ifndef CUBICLATTICECANONICALIZOR_H
#define CUBICLATTICECANONICALIZOR_H

#include <vector>
#include <iostream>
#include <array>
#include <cmath>
#include <algorithm>

extern "C" {
#include "nauty.h"
}

#include "CubicLattice.h"
#include "GraphContainer.h"
#include "VertexEmbedList.h"

class CubicLatticeCanonicalizor
{
public:
    typedef std::array<double, 3> CartesianCoords; /// TODO: encode color in the first coordinate (0 and 1 for colored vertices and 2 for uncolored vertices)
private:
     /******** private variables ********/
    GraphContainer *Container; /// pointer to container object (memory owned externally!)

    CubicLattice *Lattice; /// pointer to cubic lattice object (memory owned externally!)

    std::vector<CartesianCoords> VerticesCartesian; /// embedded vertices Cartesian coordinates on the cubic lattice

    std::vector<CartesianCoords> VerticesCartesianCOM; /// shifted

    CartesianCoords COM; /// center of mass vector

    VertexEmbedList OriginalList; /// original embedded list

    /******** private functions ********/
    void ComputeCOMandShift(); /// compute COM vector and shift coordinates of vertices to COM frame

    void AccumulateCartesianVector(CartesianCoords& a, const CartesianCoords& b); /// a += b

    void AccumulateCartesianVectorNeg(CartesianCoords& a, const CartesianCoords& b); /// a -= b

    void ScalarMultCartesianVector(CartesianCoords& a, double s); /// a *= s

    double NormSqCartesianVector(const CartesianCoords& a); /// |\vec{a}|^2

    double NormCartesianVector(const CartesianCoords &a); /// |\vec{a}|

    CartesianCoords RotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords& p, const CartesianCoords& x); /// rotate vector x about axis by amount theta, where axis is located at p

    CartesianCoords RotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords& x); /// rotate vector x about axis (through origin) by amount theta

    CartesianCoords Inversion(const CartesianCoords& p, const CartesianCoords& x); /// inversion of vector x about about p

    CartesianCoords InversionComposedWithRotationTransform(double theta, const CartesianCoords& axis, const CartesianCoords& p, const CartesianCoords& x); /// rotation followed by inversion

    std::vector<std::vector<CartesianCoords>> ApplyRotationsFourFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices); /// rotations about four-fold axis

    std::vector<std::vector<CartesianCoords>> ApplyRotationsTwoFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices); /// rotations about two-fold axis

    std::vector<std::vector<CartesianCoords>> ApplyRotationsThreeFoldAxisOnGraph(const CartesianCoords& p, const std::vector<CartesianCoords>& vertices); /// rotations about three-fold axis

    VertexEmbedList ConvertCartesianVectorToVertexEmbedList(const std::vector<CartesianCoords>& vec); /// vector of CartesianVector objects to VertexEmbedList

    void GenerateAllPermutationsWithoutRepeats(std::vector<std::vector<double>> &lists, const std::vector<double>& s, std::vector<double> &temp); /// recursive

    void GenerateAllPermutationsWithoutRepeats(std::vector<CartesianCoords> &lists, const std::vector<double>& s); /// hard-coded

    std::vector<std::vector<CartesianCoords>> PermutationsOnAllVertices(const std::vector<CartesianCoords>& vertices); /// perform permutations on all vertices of graphs

public:
    CubicLatticeCanonicalizor(GraphContainer *container, CubicLattice *lattice, const VertexEmbedList& embedList); /// constructor

    VertexEmbedList GetCanonical(); /// get the canonical graph (main interface routine)

    VertexEmbedList GetCanonicalOld(); /// old way (using rotations)

    void PrintVerticesCartesian(); /// debugging routine

    void PrintVerticesCartesianCOM(); /// debugging routine

};

/// output of list of coorindates
std::ostream& operator<<(std::ostream& os, const CubicLatticeCanonicalizor::CartesianCoords& list);

/// relational operators for CartesianCoords
bool operator<(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs);
bool operator==(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs);
bool operator!=(const CubicLatticeCanonicalizor::CartesianCoords& lhs, const CubicLatticeCanonicalizor::CartesianCoords& rhs);
bool operator<(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs);
bool operator==(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs);
bool operator!=(const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& lhs, const std::vector<CubicLatticeCanonicalizor::CartesianCoords>& rhs);
bool operator<(const std::pair<CubicLatticeCanonicalizor::CartesianCoords,int>& lhs, const std::pair<CubicLatticeCanonicalizor::CartesianCoords,int>& rhs);

#endif // CUBICLATTICECANONICALIZOR_H
