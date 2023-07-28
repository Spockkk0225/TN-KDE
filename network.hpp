#ifndef network_hpp
#define network_hpp

#include "edge.hpp"

struct Model {
    double test = 0;

    string network_filename;
    string output_filename;
    string stat_filename;

    double prepare_time, running_time;
    double dijkstra_time;

    int n, m;
    int start_time, end_time;
    double lixel_interval_length;
    vector<Edge*> edge_set;
    vector<vector<int>> network;
    vector<Lixel> lixel_set;
    unordered_map<int, unordered_map<int, int>> node1_node2_edge_index;

    // Method type & kernel function type
    Param::Method method_type;
    Param::Kernel kernel_type;
    double bandwidth;

    // Edge Mask
    Param::Prune mask_switch = Param::OFF;
//    Param::Prune real_switch = Param::OFF;
    Param::State edge_state = Param::OLD;
    vector<int> node_mask;
    vector<int> masked_node;
    vector<int> full_masked_c_edge;
    vector<int> full_masked_d_edge;
    vector<int> part_masked_edge;
    Item item1, item2, item3, item4;
    int lixel_id;
    vector<double> lixels_aggregation;

    // Dijkstra from a lixel
    Lixel cur_lixel;
    vector<SPNode> sp_lixel_vec;
    unordered_set<int> visit_edges;

    // Dijkstra from a node
    int cur_node_a;
    int cur_node_b;
    vector<SPNode> sp_node_a_vec;
    vector<SPNode> sp_node_b_vec;

    // Multiple Query Time Window
    int windows_num;
    vector<double> start_times;
    vector<double> end_times;

    int num;
    int H;
};

#endif /* network_hpp */