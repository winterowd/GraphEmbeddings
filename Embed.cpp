#include <iostream>

#include "GraphEmbedder.h"

int main(int argc, char *argv[])
{
    std::cout << "EMBEDDING GRAPHS!\n";
    GraphEmbedder MyEmbedder(argc, argv);
    auto start = std::chrono::high_resolution_clock::now();
    MyEmbedder.Embed();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> s_int = end-start;
    std::cout << s_int.count() << "s\n";
}
