#include <iostream>

#include "TriangularLattice.h"

TriangularLattice::TriangularLattice(unsigned int n) : AbstractLattice("Triangular"), N(n)
{

}

int TriangularLattice::GetSiteIndex(const std::vector<unsigned int>& indices) const
{
    return this->N*indices[1]+indices[0];
}

void TriangularLattice::GetSiteCoordinates(int index, std::vector<unsigned int>& indices)
{
    indices[0] = index%this->N;
    indices[1] = index/this->N;
}

int TriangularLattice::GetNearestNeighbor(int siteIndex, int nnIndex)
{
    std::vector<unsigned int> temp(2);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnIndex) {
    case 1:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 2:
        temp[0] = (temp[0]+this->N-1)%this->N;
        break;
    case 3:
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
     case 5:
        temp[0] = (temp[0]+1)%this->N;
        break;
     case 6:
        temp[1] = (temp[1]+1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid nearest neighbor index for TriangularLattice!");
    }
    return this->GetSiteIndex(temp);
}

int TriangularLattice::GetNextNearestNeighbor(int siteIndex, int nnIndex)
{
    std::cerr << "ERROR! GetNextNearestNeighbor not implemented for the triangular lattice!\n";
    std::exit(1);
}

int TriangularLattice::GetThirdNearestNeighbor(int siteIndex, int nnIndex)
{
    std::cerr << "ERROR! GetThirdNearestNeighbor not implemented for the triangular lattice!\n";
    std::exit(1);
}

int TriangularLattice::GetFourthNearestNeighbor(int siteIndex, int nnIndex)
{
    std::cerr << "ERROR! GetFourthNearestNeighbor not implemented for the triangular lattice!\n";
    std::exit(1);
}

bool TriangularLattice::AreNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool TriangularLattice::AreNNN(unsigned int index1, unsigned int index2)
{
    std::cerr << "ERROR! AreNNN not implemented for the triangular lattice!\n";
    std::exit(1);
}

bool TriangularLattice::AreThirdNN(unsigned int index1, unsigned int index2)
{
    std::cerr << "ERROR! AreThirdNN not implemented for the triangular lattice!\n";
    std::exit(1);
}

bool TriangularLattice::AreFourthNN(unsigned int index1, unsigned int index2)
{
    std::cerr << "ERROR! AreFourthNN not implemented for the triangular lattice!\n";
    std::exit(1);
}

int TriangularLattice::GetManhattanDistance(int neighborDepth)
{
    std::cerr << "WARNGING! TriangularLattice::GetManhattanDistance not yet implemented!\n";
    return 1;
}
