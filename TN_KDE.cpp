#include "TN_KDE.hpp"

void TN_KDE(Model &model) {
    
    for (int l = 0;l < model.lixel_set.size(); l++) {
        if (l % 1000 == 0) cout << l << " / " << model.lixel_set.size() << endl;
        
        model.cur_lixel = model.lixel_set[l];

        clock_t start = clock();
        if (model.method_type == Param::RQS) {
            Dijkstra(model, model.sp_lixel_vec);
        }   else {
            ShortestPathSharing(model, l);
        }
        clock_t end = clock();
        model.dijkstra_time += (double)(end - start) / CLOCKS_PER_SEC;
        
        GetKDE(model);
        model.lixel_set[l] = model.cur_lixel;
    }
}


void ShortestPathSharing(Model &model, int lixel_index) {
    if (lixel_index == 0 ||
        model.lixel_set[lixel_index].edge_index != model.lixel_set[lixel_index - 1].edge_index) {
        model.cur_node_a = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_1;
        model.cur_node_b = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_2;
        
        Dijkstra(model, model.sp_node_a_vec);
        Dijkstra(model, model.sp_node_b_vec);
    }
    
    for (int v = 0;v < model.n; v++) {
        model.sp_lixel_vec[v].cur_sp_value =
            min(model.sp_node_a_vec[v].cur_sp_value + model.lixel_set[lixel_index].dist_1,
                model.sp_node_b_vec[v].cur_sp_value + model.lixel_set[lixel_index].dist_2);
    }
}

void InitDijkstra(Model &model, vector<SPNode> &sp, priority_queue<SPNode, vector<SPNode>, SPNode_cmp> &pq) {
    
    for (int v = 0;v < model.n; v++) {
        sp[v].node_index = v; // ?
        sp[v].cur_sp_value = INF;
    }
    
    if (model.method_type == Param::RQS) {
        int node1 = model.edge_set[model.cur_lixel.edge_index]->node_1;
        int node2 = model.edge_set[model.cur_lixel.edge_index]->node_2;
        sp[node1].cur_sp_value = model.cur_lixel.dist_1;
        sp[node2].cur_sp_value = model.cur_lixel.dist_2;
        
        pq.push(sp[node1]);
        pq.push(sp[node2]);
        
    }   else {
        int node = model.cur_node_a;
        // search for a first, then b
        if (node == -1) node = model.cur_node_b; else model.cur_node_a = -1;
        sp[node].cur_sp_value = 0;
        pq.push(sp[node]);
    }
}

void Dijkstra(Model &model, vector<SPNode> &sp) {
    priority_queue<SPNode, vector<SPNode>, SPNode_cmp> pq;
    
    InitDijkstra(model, sp, pq);
    
    while (!pq.empty()) {
        SPNode cur_node = pq.top();
        pq.pop();
        
        // double-version judgement
        // if (sp[cur_node.node_index].cur_sp_value != cur_node.cur_sp_value)
        if (fabs(sp[cur_node.node_index].cur_sp_value - cur_node.cur_sp_value) > EPS) {
            continue;
        }
        
        if (cur_node.cur_sp_value > model.bandwidth) {
            continue;
        }
        
        for (int e = 0; e < model.network[cur_node.node_index].size(); e++) {
            int edge_index = model.network[cur_node.node_index][e];
            Edge* edge = model.edge_set[edge_index];
            
            // trick of finding the neighbor
            int new_node_index = edge->node_1 + edge->node_2 - cur_node.node_index;
            
            double new_sp_value = cur_node.cur_sp_value + edge->length;
            if (new_sp_value < sp[new_node_index].cur_sp_value) {
                sp[new_node_index].cur_sp_value = new_sp_value;
                pq.push(sp[new_node_index]);
            }
        }
    }
}



