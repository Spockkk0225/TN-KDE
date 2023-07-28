#include "TN_KDE.hpp"

void TN_KDE(Model &model) {

    for (int l = 0;l < model.lixel_set.size(); l++) {
        if (l % 1000 == 0) cout << l << " / " << model.lixel_set.size() << endl;
        if (l == 0 || model.lixel_set[l].edge_index != model.lixel_set[l - 1].edge_index) {
            model.edge_state = Param::NEW;
//            if (model.mask_switch == Param::ON && model.edge_set[model.lixel_set[l].edge_index]->length > 1.5 * model.lixel_interval_length) {
//                model.real_switch = Param::ON;
//            }   else {
//                model.real_switch = Param::OFF;
//            }
        }   else {
            model.edge_state = Param::OLD;
        }

        model.cur_lixel = model.lixel_set[l];

        clock_t start = clock();
        if (model.method_type == Param::RQS) {
            Dijkstra(model, model.sp_lixel_vec);
        }   else {
            ShortestPathSharing(model, l);
        }

        clock_t end = clock();
        model.dijkstra_time += (double)(end - start) / CLOCKS_PER_SEC;

        GetKDE(model, l);
        model.lixel_set[l] = model.cur_lixel;
    }
}

void ShortestPathMasking(Model &model, int lixel_index) {
    if (model.edge_state == Param::NEW) {
        model.cur_node_a = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_1;
        model.cur_node_b = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_2;

        model.visit_edges.clear();
        Dijkstra(model, model.sp_node_a_vec);
        Dijkstra(model, model.sp_node_b_vec);

        model.node_mask.resize(model.n, 0);
        model.full_masked_c_edge.clear();
        model.full_masked_d_edge.clear();
        model.part_masked_edge.clear();
//        int a=0, b=0, c=0, d=0;
        for (int e : model.visit_edges) {
//        for (int e = 0;e < model.edge_set.size(); e++) {
            int node_c = model.edge_set[e]->node_1;
            int node_d = model.edge_set[e]->node_2;
            double dist_a_c = model.sp_node_a_vec[node_c].cur_sp_value;
            double dist_a_d = model.sp_node_a_vec[node_d].cur_sp_value;
            double dist_b_c = model.sp_node_b_vec[node_c].cur_sp_value;
            double dist_b_d = model.sp_node_b_vec[node_d].cur_sp_value;
            double edge_len = model.edge_set[e]->length;
            double dist = min(
                    min(model.sp_node_a_vec[node_c].cur_sp_value, model.sp_node_a_vec[node_d].cur_sp_value),
                    min(model.sp_node_b_vec[node_c].cur_sp_value, model.sp_node_b_vec[node_d].cur_sp_value));
            double max_dist_q_c = min(
                    (model.sp_node_a_vec[node_c].cur_sp_value + model.sp_node_b_vec[node_c].cur_sp_value + model.edge_set[model.cur_lixel.edge_index]->length) / 2,
                     min(model.sp_node_a_vec[node_c].cur_sp_value, model.sp_node_b_vec[node_c].cur_sp_value) + model.edge_set[model.cur_lixel.edge_index]->length);
            double max_dist_q_d = min(
                    (model.sp_node_a_vec[node_d].cur_sp_value + model.sp_node_b_vec[node_d].cur_sp_value + model.edge_set[model.cur_lixel.edge_index]->length) / 2,
                    min(model.sp_node_a_vec[node_d].cur_sp_value, model.sp_node_b_vec[node_d].cur_sp_value) + model.edge_set[model.cur_lixel.edge_index]->length);

            int lixel_num = (int)((model.edge_set[model.lixel_set[lixel_index].edge_index]->length + model.lixel_interval_length / 2) / model.lixel_interval_length);
            double max_q_c_minus_q_d = -INF;
            double min_q_c_minus_q_d = INF;

            double half_c = (model.edge_set[model.lixel_set[lixel_index].edge_index]->length + dist_a_c + dist_b_c) / 2 - dist_a_c;
            double half_d = (model.edge_set[model.lixel_set[lixel_index].edge_index]->length + dist_a_d + dist_b_d) / 2 - dist_a_d;
            int middle_c = (int)((half_c + model.lixel_interval_length / 2) / model.lixel_interval_length) - 1;
            int middle_d = (int)((half_d + model.lixel_interval_length / 2) / model.lixel_interval_length) - 1;
            update(middle_c, 0, lixel_num-1, model, lixel_index, node_c, node_d, max_q_c_minus_q_d, min_q_c_minus_q_d);
            update(middle_c + 1, 0, lixel_num-1, model, lixel_index, node_c, node_d, max_q_c_minus_q_d, min_q_c_minus_q_d);
            update(middle_d, 0, lixel_num-1, model, lixel_index, node_c, node_d, max_q_c_minus_q_d, min_q_c_minus_q_d);
            update(middle_d + 1, 0, lixel_num-1, model, lixel_index, node_c, node_d, max_q_c_minus_q_d, min_q_c_minus_q_d);

            if (dist > model.bandwidth) {
//                a++;

//            }   else if (fabs(dist_a_c + edge_len - dist_a_d) < EPS && fabs(dist_b_c + edge_len - dist_b_d) < EPS && max_dist_q_c + edge_len <= model.bandwidth) {
            }   else if (max_q_c_minus_q_d - EPS <= model.edge_set[e]->min_d_o_minus_c_o && max_dist_q_c + edge_len <= model.bandwidth) {
                model.full_masked_c_edge.push_back(e);
                model.node_mask[model.edge_set[e]->node_1] = 1;
                model.node_mask[model.edge_set[e]->node_2] = 1;
//                b++;

//            }   else if (fabs(dist_a_d + edge_len - dist_a_c) < EPS && fabs(dist_b_d + edge_len - dist_b_c) < EPS && max_dist_q_d + edge_len <= model.bandwidth) {
            }   else if (min_q_c_minus_q_d + EPS >= model.edge_set[e]->max_d_o_minus_c_o && max_dist_q_d + edge_len <= model.bandwidth) {
                model.full_masked_d_edge.push_back(e);
                model.node_mask[model.edge_set[e]->node_1] = 1;
                model.node_mask[model.edge_set[e]->node_2] = 1;
//                c++;

            }   else {
                model.part_masked_edge.push_back(e);
                model.node_mask[model.edge_set[e]->node_1] = 1;
                model.node_mask[model.edge_set[e]->node_2] = 1;
//                d++;
            }
        }
//        cout << a << " " << b << " " << c << " " << d << endl;

        model.masked_node.clear();
        for (int v = 0;v < model.n; v++) {
            if (model.node_mask[v] == 1) {
                model.masked_node.push_back(v);
            }
            model.sp_lixel_vec[v].cur_sp_value = model.bandwidth + 1;
        }
    }

    for (int i = 0;i < model.masked_node.size(); i++) {
        int v = model.masked_node[i];
        model.sp_lixel_vec[v].cur_sp_value =
                min(model.sp_node_a_vec[v].cur_sp_value + model.lixel_set[lixel_index].dist_1,
                    model.sp_node_b_vec[v].cur_sp_value + model.lixel_set[lixel_index].dist_2);
    }

}

void ShortestPathSharing(Model &model, int lixel_index) {
    if (model.mask_switch == Param::ON) {
        ShortestPathMasking(model, lixel_index);
        return;
    }

    if (model.edge_state == Param::NEW) {
        model.cur_node_a = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_1;
        model.cur_node_b = model.edge_set[model.lixel_set[lixel_index].edge_index]->node_2;

        model.visit_edges.clear();
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
            model.visit_edges.insert(edge_index);

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

void update(int pos, int L, int R, Model &model, int &lixel_index, int &node_c, int &node_d, double &max_q_c_minus_q_d, double &min_q_c_minus_q_d) {
    int i = pos;
    i = max(i, L);
    i = min(i, R);
    double dist_q_c = min(model.lixel_set[lixel_index + i].dist_1 + model.sp_node_a_vec[node_c].cur_sp_value,
                              model.lixel_set[lixel_index + i].dist_2 + model.sp_node_b_vec[node_c].cur_sp_value);
    double dist_q_d = min(model.lixel_set[lixel_index + i].dist_1 + model.sp_node_a_vec[node_d].cur_sp_value,
                              model.lixel_set[lixel_index + i].dist_2 + model.sp_node_b_vec[node_d].cur_sp_value);
    max_q_c_minus_q_d = max(max_q_c_minus_q_d, dist_q_c - dist_q_d);
    min_q_c_minus_q_d = min(min_q_c_minus_q_d, dist_q_c - dist_q_d);
}

