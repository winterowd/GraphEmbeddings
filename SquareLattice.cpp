#include <iostream>

#include "SquareLattice.h"

SquareLattice::SquareLattice(unsigned int n) : AbstractLattice("Square"), N(n)
{

}

int SquareLattice::GetSiteIndex(const std::vector<unsigned int>& indices) const
{
    return this->N*indices[1]+indices[0];
}

void SquareLattice::GetSiteCoordinates(int index, std::vector<unsigned int>& indices)
{
    indices[0] = index%this->N;
    indices[1] = index/this->N;
}

int SquareLattice::GetNearestNeighbor(int siteIndex, int nnIndex)
{
    std::vector<unsigned int> temp(2);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnIndex) {
    case 1:
        temp[0] = (temp[0]+1)%this->N;
        break;
    case 2:
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 3:
        temp[0] = (temp[0]+this->N-1)%this->N;
        break;
    case 4:
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid nearest neighbor index for SquareLattice!");
    }
    return this->GetSiteIndex(temp);
}

/// all points with Euclidean distance \sqrt{2} from siteIndex
int SquareLattice::GetNextNearestNeighbor(int siteIndex, int nnnIndex)
{
    std::vector<unsigned int> temp(2);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnnIndex) {
    case 1:
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
     default:
        throw std::invalid_argument("Invalid next nearest neighbor index for SquareLattice!");
    }
    return this->GetSiteIndex(temp);
}

/// all points with Euclidean distance 2 from siteIndex
int SquareLattice::GetThirdNearestNeighbor(int siteIndex, int nnnIndex)
{
    std::vector<unsigned int> temp(2);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnnIndex) {
    case 1:
        temp[0] = (temp[0]+2)%this->N;
        break;
    case 2:
        temp[1] = (temp[1]+2)%this->N;
        break;
    case 3:
        temp[0] = (temp[0]+this->N-2)%this->N;
        break;
    case 4:
        temp[1] = (temp[1]+this->N-2)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid third nearest neighbor index for SquareLattice!");
    }
    return this->GetSiteIndex(temp);
}

/// all points with Euclidean distance \sqrt{5} from siteIndex
int SquareLattice::GetFourthNearestNeighbor(int siteIndex, int nnnIndex)
{
    std::vector<unsigned int> temp(2);
    this->GetSiteCoordinates(siteIndex, temp);
    switch (nnnIndex) {
    case 1:
        temp[0] = (temp[0]+2)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 2:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+2)%this->N;
        break;
    case 3:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+2)%this->N;
        break;
    case 4:
        temp[0] = (temp[0]+this->N-2)%this->N;
        temp[1] = (temp[1]+1)%this->N;
        break;
    case 5:
        temp[0] = (temp[0]+this->N-2)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
    case 6:
        temp[0] = (temp[0]+this->N-1)%this->N;
        temp[1] = (temp[1]+this->N-2)%this->N;
        break;
    case 7:
        temp[0] = (temp[0]+1)%this->N;
        temp[1] = (temp[1]+this->N-2)%this->N;
        break;
    case 8:
        temp[0] = (temp[0]+2)%this->N;
        temp[1] = (temp[1]+this->N-1)%this->N;
        break;
     default:
        throw std::invalid_argument("Invalid fourth nearest neighbor index for SquareLattice!");
    }
    return this->GetSiteIndex(temp);
}

bool SquareLattice::AreNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool SquareLattice::AreNNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetNextNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool SquareLattice::AreThirdNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetThirdNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}

bool SquareLattice::AreFourthNN(unsigned int index1, unsigned int index2)
{
    for (int i=1; i<=this->NbrNN; ++i)
        if (this->GetFourthNearestNeighbor(index1, i) == index2)
            return true;
    return false;
}
