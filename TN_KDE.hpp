#ifndef TN_KDE_hpp
#define TN_KDE_hpp

#include "kernel.hpp"
#include "library.hpp"

void TN_KDE(Model &model);

void ShortestPathSharing(Model &model, int lixel_index);
void Dijkstra(Model &model, vector<SPNode> &sp);
void InitDijkstra(Model &model, vector<SPNode> &sp, priority_queue<SPNode, vector<SPNode>, SPNode_cmp> &pq);
void update(int pos, int L, int R, Model &model, int &lixel_index, int &node_c, int &node_d, double &max_q_c_minus_q_d, double &min_q_c_minus_q_d);

#endif /* TN_KDE_hpp */