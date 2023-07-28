#include "init.hpp"

Param::Method CheckMethod(char* arg) {
    if (strcmp(arg, "RQS") == 0) { return Param::RQS; }
    if (strcmp(arg, "SPS") == 0) { return Param::SPS; }
    if (strcmp(arg, "RAS") == 0) { return Param::RAS; }
    if (strcmp(arg, "RTS") == 0) { return Param::RTS; }
    if (strcmp(arg, "RFS") == 0) { return Param::RFS; }
    if (strcmp(arg, "DRFS") == 0) { return Param::DRFS; }
    if (strcmp(arg, "ADA") == 0) { return Param::ADA; }

    return Param::RQS;
}

Param::Kernel CheckType(char* arg) {
    if (strcmp(arg, "Triangular") == 0) { return Param::Triangular; }
    return Param::Triangular;
}

Param::Prune CheckSwitch(char* arg) {
    if (strcmp(arg, "ON") == 0) { return Param::ON; }
    return Param::OFF;
}

void InitParameters(Model &model, int argc, char** argv) {
    model.network_filename = argv[1];
    model.output_filename = argv[2];
    model.stat_filename = argv[3];
    model.method_type = CheckMethod(argv[4]);
    model.kernel_type = CheckType(argv[5]);
    model.mask_switch = CheckSwitch(argv[6]);
    model.lixel_interval_length = atoi(argv[7]);
    model.bandwidth = atoi(argv[8]);
    model.H = atoi(argv[9]);
    if (argc > 10) {
        model.windows_num = atoi(argv[10]);
        for (int i = 0;i < model.windows_num; i++) {
            model.start_times.push_back(atof(argv[i * 2 + 11]));
            model.end_times.push_back(atof(argv[i * 2 + 11 + 1]));
        }
    } else {
        model.windows_num = 0;
        model.start_time = -1;
        model.end_time = -1;
    }

//    model.network_filename = "network.SanFrancisco.txt";
//    model.output_filename = "RFS.SanFrancisco.txt";
//    model.stat_filename = "statistic.txt";
//    model.method_type = Param::RFS;
//    model.kernel_type = Param::Triangular;
//    model.mask_switch = Param::ON;
//    model.lixel_interval_length = 50;
//    model.bandwidth = 3000;
//    model.H = 0;
//    model.windows_num = 6;
//    model.start_times = vector<double>{0.11, 0.12, 0.13, 0.14, 0.15, 0.16};
//    model.end_times = vector<double>{0.81, 0.82, 0.83, 0.84, 0.85, 0.86};

    cout << "init_parameters finished" << endl;
}

Edge* NewEdge(Model &model) {
    if (model.method_type == Param::RQS || model.method_type == Param::SPS)  return new Baseline;
    if (model.method_type == Param::RAS) return new RAS;
    if (model.method_type == Param::RTS) return new RTS;
    if (model.method_type == Param::RFS) return new RFS;
    if (model.method_type == Param::DRFS) return new DRFS;
    if (model.method_type == Param::ADA) return new ADA;
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
    }

    sort(p.begin(), p.end(), PointDistCmp);

    for (int i = 0;i < p.size(); i++) {
        model.edge_set[p[i].edge_index]->min_d_o_minus_c_o = min(model.edge_set[p[i].edge_index]->min_d_o_minus_c_o,
                                                                 p[i].dist_2 - p[i].dist_1);
        model.edge_set[p[i].edge_index]->max_d_o_minus_c_o = max(model.edge_set[p[i].edge_index]->max_d_o_minus_c_o,
                                                                 p[i].dist_2 - p[i].dist_1);
        model.edge_set[p[i].edge_index]->point_set.push_back(p[i]);
    }

    cout << "load_network finished : event = " << point_num << endl;

    // init for Dijkstra array
    model.sp_lixel_vec.resize(model.n, SPNode{});
    model.sp_node_a_vec.resize(model.n, SPNode{});
    model.sp_node_b_vec.resize(model.n, SPNode{});

    if (model.windows_num == 0) { // default
        model.windows_num = 1;
        double interval = max_time - min_time;
        model.start_times.push_back(min_time + 0.2 * interval);
        model.end_times.push_back(max_time - 0.2 * interval);
        cout << model.start_times[0] << " " << model.end_times[0] << endl;
    }   else {
        sort(p.begin(), p.end(), PointTimeCmp);
        for (int i = 0;i < model.start_times.size(); i++) {
            model.start_times[i] = p[(int)(model.start_times[i] * (int)p.size())].time;
            model.end_times[i] = p[(int)(model.end_times[i] * (int)p.size())].time;
            cout << model.start_times[i] << " ==> " << model.end_times[i] << endl;
        }
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
    double cur_dist = model.lixel_interval_length / 2;
    double length = model.edge_set[edge_index]->length;
    Lixel lixel;

    while (cur_dist < length) {
        lixel.dist_1 = cur_dist;
        lixel.dist_2 = length - cur_dist;
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
    if (method_type == Param::ADA) return "ADA";
    return "Unknown";
}

string OutputKernel(Param::Kernel kernel_type) {
    if (kernel_type == Param::Triangular) return "Triangular";
    if (kernel_type == Param::Cosine) return "Cosine";
    if (kernel_type == Param::Exponential) return "Exponential";
    return "Unknown";
}

string OutputPrune(Param::Prune prune_type) {
    if (prune_type == Param::ON) return "ON";
    if (prune_type == Param::OFF) return "OFF";
    return "Unknown";
}

void Statistic(Model &model) {
    ofstream stat_file;
    stat_file.open(model.stat_filename, ios::app);

    stat_file << "Dataset=" << model.network_filename << endl;
    stat_file << "MethodType=" << OutputMethod(model.method_type) << endl;
    stat_file << "KernelType=" << OutputKernel(model.kernel_type) << endl;
    stat_file << "LixelInterval=" << model.lixel_interval_length << endl;
    stat_file << "BandWidth=" << model.bandwidth << endl;

    stat_file << "PruneType=" << OutputPrune(model.mask_switch) << endl;
    stat_file << "TimeWindow=" << model.windows_num << endl;
    stat_file << "RunningTime=" << model.running_time << endl;
    stat_file << endl;

    stat_file.close();
}

void Visualize(Model &model, int id) {
    ofstream output_file;
    string ch = ".?";
    ch[1] = id + '0';
    output_file.open(model.output_filename + ch);

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

void InitEdgeStructure(Model &model, int i) {
    model.start_time = model.start_times[i];
    model.end_time = model.end_times[i];
    double b_s = model.bandwidth;
    double b_t = (model.end_time - model.start_time) / 2.0;
    for (Edge* e : model.edge_set) {
        e->b_s = b_s;
        e->b_t = b_t;
    }

    if (model.method_type == Param::RQS || model.method_type == Param::SPS) return;

    for (int e = 0;e < model.m; e++) {
        if (model.method_type == Param::RFS || model.method_type == Param::DRFS || model.method_type == Param::ADA) {
            model.edge_set[e]->SetTime(model.start_time, model.end_time);
        }
        if (i == 0) model.edge_set[e]->Init(model.kernel_type); // only initialized for the first time
        model.edge_set[e]->SubInit();
        if (model.method_type == Param::DRFS) {
            model.edge_set[e]->SetH(model.H);
        }
    }
}
