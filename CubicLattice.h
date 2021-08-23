#ifndef CUBICLATTICE_H
#define CUBICLATTICE_H

#include "AbstractLattice.h"

class CubicLattice : public AbstractLattice
{
private:
    unsigned int N; /// dimension of lattice

    unsigned int NSq;

public:
    CubicLattice(unsigned int n);

    ~CubicLattice() {};

    int GetSiteIndex(const std::vector<unsigned int>& indices) const override;

    void GetSiteCoordinates(int index, std::vector<unsigned int>& indices) override;

    int GetNearestNeighbor(int siteIndex, int nnIndex) override;

    int GetNextNearestNeighbor(int siteIndex, int nnIndex) override;

    unsigned int GetN() const override { return this->N; }

    unsigned int GetNbrNN() const override { return this->NbrNN; }

    unsigned int GetNbrNNN() const override { return this->NbrNNN; }

    unsigned int GetDim() const override { return this->Dim; }

    bool AreNN(unsigned int index1, unsigned int index2) override;

    bool AreNNN(unsigned int index1, unsigned int index2) override;

    static const unsigned int Dim = 3;

    static const unsigned int NbrNN = 6; /// vectors labeled by 1,2,...,NbrNN

    static const unsigned int NbrNNN = 12; /// vectors labeled by 1,2,...,NbrNNN

    int GetOppositeNNDirection(int nnIndex) { return (nnIndex+2)%this->NbrNN+1; }
};

#endif // CUBICLATTICE_H
