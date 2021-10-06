#include <iostream>

#include "GraphEmbedder.h"

int main(int argc, char *argv[])
{
    std::cout << "EMBEDDING GRAPHS!\n";
    GraphEmbedder MyEmbedder(argc, argv);
    auto start = std::chrono::high_resolution_clock::now();
    MyEmbedder.Embed();
    //MyEmbedder.TestInitialRootedGraphList("test_fifth_order_dense_g6.dat");
    //MyEmbedder.TestInitialRootedGraphList("test_fifth_order_dense_g6.dat" /* fixed by hand containing polymer */, 1 /* v1 */, 5 /* v2 */); /// run with fifth order, max embedding NNN, and correlator distance to NN
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> s_int = end-start;
    std::cout << s_int.count() << "s\n";
}
