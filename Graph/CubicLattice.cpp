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


/// all points with Euclidean distance \sqrt{2} from siteIndex
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


/// all points with Euclidean distance 2 from siteIndex
int CubicLattice::GetThirdNearestNeighbor(int siteIndex, int n3Index)
{
    std::vector<unsigned int> temp(3);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (n3Index) {
    case 1:
        temp[0] = (temp[0]+2)%this->N;
        break;
    case 2:
        temp[1] = (temp[1]+2)%this->N;
        break;
    case 3:
        temp[2] = (temp[2]+2)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+this->N-2)%this->N;
        break;
    case 5:
        temp[1] = (temp[1]+this->N-2)%this->N;
        break;
    case 6:
        temp[2] = (temp[2]+this->N-2)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid third nearest neighbor index for CubicLattice!");
    }
    return this->GetSiteIndex(temp);
}

/// all points with Euclidean distance \sqrt{3} from siteIndex
int CubicLattice::GetFourthNearestNeighbor(int siteIndex, int n4Index)
{
    std::vector<unsigned int> temp(3);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (n4Index) {
    case 1:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 2:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 3:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+1)%this->N;
        break;
    case 5:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 6:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 7:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
    case 8:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        temp[2] = (temp[2]+this->N-1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid fourth nearest neighbor index for CubicLattice!");
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
    for (int i=1; i<=this->NbrNNN; ++i)
        if (this->GetNextNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool CubicLattice::AreThirdNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrThirdNN; ++i)
        if (this->GetThirdNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool CubicLattice::AreFourthNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrFourthNN; ++i)
        if (this->GetFourthNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

int CubicLattice::GetManhattanDistance(int neighborDepth)
{
    int result;
    switch (neighborDepth) {
    case 1: /// nearest neighbor
        result=1;
        break;
    case 2: /// next nearest neighbor
        result=2;
        break;
    case 3: /// third nearest neighbor
        result=2;
        break;
    case 4: /// fourth nearest neighbor
        result=3;
        break;
    default:
        throw std::invalid_argument("Invalid neighborDepth for CubicLattice::GetManhattanDistance!");
    }
    return result;
}

