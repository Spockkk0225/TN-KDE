#ifndef init_hpp
#define init_hpp

#include "network.hpp"
#include "library.hpp"

void InitParameters(Model &model, int argc, char** argv);
void LoadNetwork(Model &model);
void InitEdgeStructure(Model &model, int id);
void AddLixels(Model &model);
void AddLixel(int edge_index, Model &model);
void Visualize(Model &model, int id);
void Statistic(Model &model);

Edge* NewEdge(Model &model);

#endif /* init_hpp */