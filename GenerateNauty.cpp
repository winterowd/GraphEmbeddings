#include <iostream>

#include "GraphGeneratorNauty.h"

int main(int argc, char *argv[])
{
    std::cout << "GENERATING GRAPHS USING NAUTY!\n";
    GraphGeneratorNauty MyGenerator(argc, argv);
    MyGenerator.Generate();
    std::cout << "GENERATION COMPLETE!\n";
}
