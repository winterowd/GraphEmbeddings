# GraphEmbeddings
Generation of connected graphs via NAUTY and their embeddings on Bravais lattices. 

g++ -o GenerateNauty gtools.o geng.o -O3 -std=c++11 -DDEBUG -DGENG_MAIN=geng_main -I/Users/christopherwinterowd/sandbox/graph_embeddings/nauty27r1 GenerateNauty.cpp GenericConnectedGraph.cpp GraphGenerator.cpp GraphGeneratorNauty.cpp GraphContainer.cpp GraphEmbedder.cpp SquareLattice.cpp CubicLattice.cpp TriangularLattice.cpp /Users/christopherwinterowd/sandbox/graph_embeddings/nauty27r1/nauty.a -lboost_program_options

g++ -o Embed gtools.o geng.o -O3 -std=c++11 -DGENG_MAIN=geng_main -I/Users/christopherwinterowd/sandbox/graph_embeddings/nauty27r1 Embed.cpp GenericConnectedGraph.cpp GraphGenerator.cpp GraphGeneratorNauty.cpp GraphContainer.cpp GraphEmbedder.cpp SquareLattice.cpp CubicLattice.cpp TriangularLattice.cpp /Users/christopherwinterowd/sandbox/graph_embeddings/nauty27r1/nauty.a -lboost_program_options
