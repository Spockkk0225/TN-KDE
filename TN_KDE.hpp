#ifndef TN_KDE_hpp
#define TN_KDE_hpp

#include "kernel.hpp"

void TN_KDE(Model &model);

void ShortestPathSharing(Model &model, int lixel_index);
void Dijkstra(Model &model, vector<SPNode> &sp);
void InitDijkstra(Model &model, vector<SPNode> &sp, priority_queue<SPNode, vector<SPNode>, SPNode_cmp> &pq);


#endif /* TN_KDE_hpp */
