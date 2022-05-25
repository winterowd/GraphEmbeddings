#include <iostream>

#include "CanonicalGraphManager.h"
#include "WrapperRoutinesForEmbedding.h"
#include "XiRecursionRooted.h"
#include "CorrelatorParameters.h"

int main(int argc, char *argv[])
{
    CorrelatorParameters MyParameters(argc, argv);

    int maxManhattanDistance = MyParameters.GetMaxManhattanDistance();
    int maxBonds = maxManhattanDistance;
    CanonicalGraphManager MyManager(maxBonds);

    CubicLattice MyLattice(100);

    GiNaC::ex resultExpanded, resultFull;
    for (int l=1; l<=maxBonds; ++l)
    {
        std::cout << "Processing graphs with " << l << " bonds!\n";
        for (int i=0; i<MyManager.GetNbrRootedGraphs(l,2); ++i)
        {
            auto tempContainer = MyManager.GetRootedGraph(l,i,2);
            std::cout << "Processing graph number " << i+1 << "\n";
            /*std::cout << "GRAPH: (" << l << "," << i << ")\n";
            std::cout << tempContainer;*/
            auto resultEmbedding = WrapperRouintesForEmbedding::ComputeRootedCanonicalEmbeddingsAndCountsCubic(tempContainer, MyParameters.GetMaxEmbeddingLength(), MyParameters.GetCorrelatorLength(), maxManhattanDistance).second;
            auto numBondCombos = resultEmbedding.size();
            for (int j=0; j<numBondCombos; ++j)
            {
                auto currentCanonReps = std::get<0>(resultEmbedding[j]);
                auto currentCanonCounts = std::get<1>(resultEmbedding[j]);
                auto currentBondCombo = std::get<2>(resultEmbedding[j]);
                /*std::cout << "BOND_COMBO: (";
                for (int k=0; k<currentBondCombo.size(); ++k)
                    std::cout << currentBondCombo[k] << ",";
                std::cout << ")\n";*/
#ifdef DEBUG
                if (currentCanonReps.size()!=currentCanonCounts.size())
                    throw std::invalid_argument("VERTEX EMBED LISTS DIFFERENT FROM COUNTS!\n");
#endif
                for (int k=0; k<currentCanonCounts.size(); ++k)
                {
                    /*std::cout << "EMBED_CANON_" << k << ":\n";
                    std::cout << currentCanonReps[k];
                    std::cout << "EMBEDDING: " << currentCanonCounts[k] << "\n";*/
                    if (MyParameters.IncludeStaticQuarks()) /// TODO: multiple flavors
                    {
                        XiRecursionRooted<ZClusterStaticQuarkArbEmbedding> tempXi(&MyManager, tempContainer, currentCanonReps[k], &MyLattice, currentCanonCounts[k], maxManhattanDistance, MyParameters.GetMaxOrderH1(), MyParameters.GetMaxOrderHBar1());
                        resultExpanded += tempXi.GetExpandedXiGiNaCWithCoefficient();
                    }
                    else
                    {
                        XiRecursionRooted<ZClusterPureGaugeArbEmbedding> tempXi(&MyManager, tempContainer, currentCanonReps[k], &MyLattice, currentCanonCounts[k], maxManhattanDistance);
                        resultExpanded += tempXi.GetExpandedXiGiNaCWithCoefficient();
                    }
                    //std::cout << "DEBUG_EXPANDED: " << resultExpanded << "\n";
                    //resultFull += tempXi.GetFullXiGiNaCWithCoefficient();
                }
            }
        }
        std::cout << "EXPANDED_RESULT_" << l << ": " << resultExpanded << "\n";
        //std::cout << "FULL_RESULT_" << l << ": " << resultFull << "\n";
    }
    std::cout << "FINAL_EXPANDED_RESULT: " << resultExpanded << "\n";
    //std::cout << "FINAL_FULL_RESULT: " << resultFull << "\n";
    return 0;
}
