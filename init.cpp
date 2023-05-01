#include "init.hpp"

Param::Method CheckMethod(char* arg) {
    if (strcmp(arg, "RQS") == 0) { return Param::RQS; }
    if (strcmp(arg, "SPS") == 0) { return Param::SPS; }
    if (strcmp(arg, "RAS") == 0) { return Param::RAS; }
    if (strcmp(arg, "RTS") == 0) { return Param::RTS; }
    if (strcmp(arg, "RFS") == 0) { return Param::RFS; }
    if (strcmp(arg, "DRFS") == 0) { return Param::DRFS; }

    return Param::RQS;
}

Param::Kernel CheckType(char* arg) {
    if (strcmp(arg, "Triangular") == 0) { return Param::Triangular; }
    return Param::Triangular;
}

void InitParameters(Model &model, int argc, char** argv) {
   model.network_filename = argv[1];
   model.output_filename = argv[2];
   model.stat_filename = argv[3];
   model.method_type = CheckMethod(argv[4]);
   model.kernel_type = CheckType(argv[5]);
   model.lixel_interval_length = atoi(argv[6]);
   model.bandwidth = atoi(argv[7]);
   model.H = atoi(argv[8]);
   if (argc > 9) {
       model.start_time = atoi(argv[9]);
       model.end_time = atoi(argv[10]);
   } else {
       model.start_time = -1;
       model.end_time = -1;
   }

//    model.network_filename = "network.JohnsCreek.txt";
//    model.output_filename = "output.JohnsCreek.txt";
//    model.stat_filename = "statistic.txt";
//    model.method_type = Param::RFS;
//    model.kernel_type = Param::Triangular;
//    model.lixel_interval_length = 50;
//    model.bandwidth = 1000;
//    model.H = 10;
//    model.start_time = -1;
//    model.end_time = -1;

    cout << "init_parameters finished" << endl;
}

Edge* NewEdge(Model &model) {
    if (model.method_type == Param::RQS || model.method_type == Param::SPS)  return new Baseline;
    if (model.method_type == Param::RAS) return new RAS;
    if (model.method_type == Param::RTS) return new RTS;
    if (model.method_type == Param::RFS) return new RFS;
    if (model.method_type == Param::DRFS) return new DRFS;
    return nullptr;
}


void LoadNetwork(Model &model) {
    ifstream network_file;
    network_file.open(model.network_filename);
    
    network_file >> model.n;
    network_file >> model.m;
    model.network.resize(model.n, vector<int>(0));
   // model.edge_set.resize(model.m);
    
    int u, v;
    double l1, l2, t;
    for (int i = 0;i < model.m; i++) {
        Edge* e = NewEdge(model);
        network_file >> e->node_1;
        network_file >> e->node_2;
        network_file >> e->length;
        e->node_1--;
        e->node_2--;
        model.edge_set.push_back(e);
        model.node1_node2_edge_index[e->node_1][e->node_2] = i;
        model.network[e->node_1].push_back(i);
        model.network[e->node_2].push_back(i);
    }
    
    vector<Point> p;
    Point temp_p;
    int point_num = 0;
    double min_time = 1000000000, max_time = 0;
    while (network_file >> u >> v) {
        u--; v--;
        network_file >> l1 >> l2 >> t;
        temp_p.edge_index = model.node1_node2_edge_index[u][v];
        temp_p.dist_1 = l1;
        temp_p.dist_2 = l2;
        temp_p.time = t / 100;
        
        // there are some wrong record in NewYork dataset
        if (model.network_filename == "network.NewYork.txt" && temp_p.time > 14515776) continue;

        min_time = min(min_time, temp_p.time);
        max_time = max(max_time, temp_p.time);
        p.push_back(temp_p);
        point_num++;
        
        // shift to generate more events
        // if (temp_p.dist_2 > 1000000000) {
        //     for (int i = 0;i < 60; i++) {
        //         temp_p.dist_1 += 0.1;
        //         temp_p.dist_2 -= 0.1;
        //         temp_p.time += 10000;
        //         p.push_back(temp_p);
        //     }
        // }
    }
    
    sort(p.begin(), p.end(), PointDistCmp);
    
    for (int i = 0;i < p.size(); i++) {
        model.edge_set[p[i].edge_index]->point_set.push_back(p[i]);
    }
    
    cout << "load_network finished : event = " << point_num << endl;
    
    // init for Dijkstra array
    model.sp_lixel_vec.resize(model.n, SPNode{});
    model.sp_node_a_vec.resize(model.n, SPNode{});
    model.sp_node_b_vec.resize(model.n, SPNode{});
    
    if (model.start_time < 0 && model.end_time < 0) { // default
        double interval = max_time - min_time;
        model.start_time = min_time + 0.2 * interval;
        model.end_time = max_time - 0.2 * interval;
    }
    
    double b_s = model.bandwidth;
    double b_t = (model.end_time - model.start_time) / 2.0;
    for (Edge* e : model.edge_set) {
        e->b_s = b_s;
        e->b_t = b_t;
    }
    
    network_file.close();
}

void AddLixels(Model &model) {
    for (int e = 0;e < model.m; e++) {
        AddLixel(e, model);
    }
    
    cout << "add_lixels finished" << endl;
}

void AddLixel(int edge_index, Model &model) {
    double cur_dist = 0;
    double middle_dist;
    double next_dist;
    double length = model.edge_set[edge_index]->length;
    Lixel lixel;
    
    while (cur_dist < length) {
        next_dist = cur_dist + model.lixel_interval_length;
        if (next_dist > length) {
            next_dist = length;
        }
        
        middle_dist = (cur_dist + next_dist) / 2;
        lixel.dist_1 = middle_dist;
        lixel.dist_2 = length - middle_dist;
        lixel.edge_index = edge_index;
        model.lixel_set.push_back(lixel);
        
        cur_dist += model.lixel_interval_length;
        
     //   break;
    }
}

string OutputMethod(Param::Method method_type) {
    if (method_type == Param::RQS) return "RQS";
    if (method_type == Param::SPS) return "SPS";
    if (method_type == Param::RAS) return "RAS";
    if (method_type == Param::RTS) return "RTS";
    if (method_type == Param::RFS) return "RFS";
    return "Unknown";
}

string OutputKernel(Param::Kernel kernel_type) {
    if (kernel_type == Param::Triangular) return "Triangular";
    if (kernel_type == Param::Cosine) return "Cosine";
    if (kernel_type == Param::Exponential) return "Exponential";
    return "Unknown";
}

void Visualize(Model &model) {
    ofstream stat_file;
    stat_file.open(model.stat_filename, ios::app);
    
    stat_file << "Dataset=" << model.network_filename << endl;
    stat_file << "MethodType=" << OutputMethod(model.method_type) << endl;
    stat_file << "KernelType=" << OutputKernel(model.kernel_type) << endl;
    stat_file << "LixelInterval=" << model.lixel_interval_length << endl;
    stat_file << "BandWidth=" << model.bandwidth << endl;
    
    stat_file << "PrepareTime=" << model.prepare_time << endl;
    stat_file << "RunningTime=" << model.running_time << endl;
    stat_file << "DijkstraTime=" << model.dijkstra_time << endl;
    stat_file << endl;
    
    stat_file.close();
    
    ofstream output_file;
    output_file.open(model.output_filename);
    
    double minV = 10000000, maxV = 0;
    for (int l = 0;l < model.lixel_set.size(); l++) {
        minV = min(minV, model.lixel_set[l].KDE_value);
        maxV = max(maxV, model.lixel_set[l].KDE_value);
    }
    
    output_file << model.lixel_set.size() << endl;
    for (int l = 0;l < model.lixel_set.size(); l++) {
        model.lixel_set[l].KDE_value = (model.lixel_set[l].KDE_value - minV) / (maxV - minV);
        output_file << l << " "
                    << model.lixel_set[l].edge_index << " "
                    << model.lixel_set[l].dist_1 << " "
                    << model.lixel_set[l].dist_2 << " "
                    << model.lixel_set[l].KDE_value << endl;
    }
    
    output_file.close();
}

void InitEdgeStructure(Model &model) {
    if (model.method_type == Param::RQS || model.method_type == Param::SPS) return;

    for (int e = 0;e < model.m; e++) {
        if (model.method_type == Param::RFS || model.method_type == Param::DRFS) {
            model.edge_set[e]->SetTime(model.start_time, model.end_time);
        }
        model.edge_set[e]->Init(model.kernel_type);
        if (model.method_type == Param::DRFS) {
            model.edge_set[e]->SetH(model.H);
        }
    }
}
