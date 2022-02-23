#include <iostream>

#include "ZClusterPureGauge.h"
#include "ZClusterPureGaugeArbEmbedding.h"
#include "CanonicalGraphManager.h"
#include "TwoPointCorrelator.h"
#include "WrapperRoutinesForEmbedding.h"

int main(int argc, char *argv[])
{
    int lMax = 3;
    CanonicalGraphManager MyManager(lMax);

    CubicLattice MyLattice(100);

    for (int l=1; l<=lMax; ++l)
    {
        for (int i=0; i<MyManager.GetNbrRootedGraphs(l,2); ++i)
        {
            auto tempContainer = MyManager.GetRootedGraph(l,i,2);
            auto resultEmbedding = WrapperRouintesForEmbedding::ComputeRootedCanonicalEmbeddingsAndCountsCubicNN(tempContainer, MaxInteractionLength::NearestNeighbor);
#ifdef DEBUG
            if (std::get<1>(resultEmbedding).size()!=std::get<2>(resultEmbedding).size())
                throw std::invalid_argument("VERTEX EMBED LISTS DIFFERENT FROM COUNTS!\n");
#endif
            std::cout << "GRAPH (" << l << "," << i << "):\n";
            std::cout << tempContainer;
            auto embedLists = std::get<1>(resultEmbedding);
            auto counts = std::get<2>(resultEmbedding);
            for (int j=0; j<counts.size(); ++j)
            {

                std::cout << "NONZERO_EMBEDDING: " << j << " with counts " << counts[j] << "\n";
                std::cout << embedLists[j];
                TwoPointCorrelator tempCorrelator(tempContainer, embedLists[j], &MyLattice);
                std::cout << "CORRELATOR!\n";
                tempCorrelator.PrintCorrelatorTerms();
            }
        }
    }

    return 0;
}
