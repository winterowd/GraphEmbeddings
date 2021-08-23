#include "CubicLattice.h"

CubicLattice::CubicLattice(unsigned int n) : AbstractLattice("Cubic"), N(n), NSq(n*n)
{

}

int CubicLattice::GetSiteIndex(const std::vector<unsigned int>& indices) const
{
    return this->NSq*indices[2]+this->N*indices[1]+indices[0];
}

void CubicLattice::GetSiteCoordinates(int index, std::vector<unsigned int>& indices)
{
    indices[2] = index/this->NSq;
    auto temp = index%this->NSq;
    indices[1] = temp/this->N;
    indices[0] = temp%this->N;
}

int CubicLattice::GetNearestNeighbor(int siteIndex, int nnIndex)
{
    std::vector<unsigned int> temp(3);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnIndex) {
    case 1:
        temp[0] = (temp[0]+1)%this->N;
        break;
    case 2:
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 3:
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+this->N-1)%this->N;
        break;
    case 5:
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
    case 6:
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid nearest neighbor index for CubicLattice!");
    }
    return this->GetSiteIndex(temp);
}

int CubicLattice::GetNextNearestNeighbor(int siteIndex, int nnIndex)
{
    std::vector<unsigned int> temp(3);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnIndex) {
    case 1: /// first four are in x-y plane
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 2:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 3:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
    case 5: /// second four have one step in +z direction
        temp[0] = (temp[0]+1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 6:
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 7:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 8:
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 9: /// last four have one step in -z-direction
        temp[0] = (temp[0]+1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 10:
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 11:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 12:
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid next nearest neighbor index for CubicLattice!");
    }
    return this->GetSiteIndex(temp);
}

bool CubicLattice::AreNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool CubicLattice::AreNNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetNextNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}


