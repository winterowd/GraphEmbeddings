#ifndef GRAPHGENERATORNAUTY_H
#define GRAPHGENERATORNAUTY_H

class GraphGeneratorParametersNauty /// graph generator parameters (get from user using Boost program options)
{
private:
    unsigned int N; /// max order of vertices

    bool Disconnected; /// do we allow disconnected diagrams

    bool ProcessCommandLine(int argc, char *argv[]);

public:

    GraphGeneratorParametersNauty(int argc, char *argv[]);

    unsigned int GetN() const { return this->N; }

    bool AllowDisconnected() const { return this->Disconnected; }

};

class GraphGeneratorNauty
{
private:
    GraphGeneratorParametersNauty Parameters;

public:
    GraphGeneratorNauty(int argc, char *argv[]);

    void Generate();
};

#endif // GRAPHGENERATORNAUTY_H
