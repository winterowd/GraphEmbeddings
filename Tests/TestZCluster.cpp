#include <iostream>

#include "ZClusterPureGauge.h"
#include "ZClusterPureGaugeArbEmbedding.h"
#include "CanonicalGraphManager.h"
#include "TwoPointCorrelator.h"
#include "WrapperRoutinesForEmbedding.h"

int main(int argc, char *argv[])
{
    int lMax = 4;
    CanonicalGraphManager MyManager(lMax);

    CubicLattice MyLattice(100);

    for (int l=1; l<=lMax; ++l)
    {
        for (int i=0; i<MyManager.GetNbrRootedGraphs(l,2); ++i)
        {
            auto tempContainer = MyManager.GetRootedGraph(l,i,2);
            auto resultEmbedding = WrapperRouintesForEmbedding::ComputeRootedCanonicalEmbeddingsAndCountsCubic(tempContainer, MaxInteractionLength::NearestNeighbor, MaxInteractionLength::NearestNeighbor, l).second;
            auto numBondCombos = resultEmbedding.size();
            for (int j=0; j<numBondCombos; ++j)
            {
                auto currentCanonReps = std::get<0>(resultEmbedding[j]);
                auto currentCanonCounts = std::get<1>(resultEmbedding[j]);
                auto currentBondCombo = std::get<2>(resultEmbedding[j]);
                std::cout << "BOND_COMBO: (";
                for (int k=0; k<currentBondCombo.size(); ++k)
                    std::cout << currentBondCombo[k] << ",";
                std::cout << ")\n";
#ifdef DEBUG
                if (currentCanonReps.size()!=currentCanonCounts.size())
                    throw std::invalid_argument("VERTEX EMBED LISTS DIFFERENT FROM COUNTS!\n");
#endif
                std::cout << "GRAPH (" << l << "," << i << "):\n";
                std::cout << tempContainer;
                std::cout << "NUMBER_OF_CANONICAL_EMBEDDINGS: " << currentCanonCounts.size() << "\n";
                for (int k=0; k<currentCanonCounts.size(); ++k)
                {

                    std::cout << "NONZERO_EMBEDDING: " << k << " with counts " << currentCanonCounts[k] << "\n";
                    std::cout << currentCanonReps[k];
                    TwoPointCorrelator<ZClusterPureGaugeArbEmbedding> tempCorrelator(tempContainer, currentCanonReps[k], &MyLattice);
                    std::cout << "CORRELATOR!\n";
                    tempCorrelator.PrintCorrelatorTerms();
                    std::cout << "GiNaC_FULLCORRELATOR!\n";
                    std::cout << tempCorrelator.GetFullCorrelatorGiNaC() << "\n";
                    std::cout << "GiNaC_EXPANDEDCORRELATOR!\n";
                    std::cout << tempCorrelator.GetExpandedCorrelatorGiNaC() << "\n";
                }
            }
        }
    }

    return 0;
}
