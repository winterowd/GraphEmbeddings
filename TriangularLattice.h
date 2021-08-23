#ifndef TRIANGULARLATTICE_H
#define TRIANGULARLATTICE_H

#include "AbstractLattice.h"

class TriangularLattice : public AbstractLattice
{
private:
    unsigned int N; /// dimension of lattice

public:
    TriangularLattice(unsigned int n);

    ~TriangularLattice() {};

    int GetSiteIndex(const std::vector<unsigned int>& indices) const override;

    void GetSiteCoordinates(int index, std::vector<unsigned int>& indices) override;

    int GetNearestNeighbor(int siteIndex, int nnIndex) override;

    int GetNextNearestNeighbor(int siteIndex, int nnnIndex) override;

    unsigned int GetN() const override { return this->N; }

    unsigned int GetNbrNN() const override { return this->NbrNN; }

     unsigned int GetNbrNNN() const override { return this->NbrNNN; }

    unsigned int GetDim() const override { return this->Dim; }

    bool AreNN(unsigned int index1, unsigned int index2) override;

    bool AreNNN(unsigned int index1, unsigned int index2) override;

    static const unsigned int Dim = 2;

    static const unsigned int NbrNN = 6;  /// vectors labeled by 1,2,...,NbrNN

    static const unsigned int NbrNNN = 4;  /// vectors labeled by 1,2,...,NbrNNN

    int GetOppositeNNDirection(int nnIndex) { return (nnIndex+1)%this->NbrNN+1; }
};

#endif // TRIANGULARLATTICE_H
