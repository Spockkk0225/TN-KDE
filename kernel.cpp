#include "kernel.hpp"

void GetKDE(Model &model, int lixel_index) {
    double edge_KDE_value = 0;

    if (model.mask_switch == Param::ON) {

        if (model.edge_state == Param::NEW) {
            model.lixel_id = -1;
            InitItemMasking(model);
            model.edge_state = Param::OLD;
        }

        model.lixel_id++;
        edge_KDE_value = model.lixels_aggregation[model.lixel_id];

        for (int i = 0;i < model.part_masked_edge.size(); i++) {
            int e = model.part_masked_edge[i];
            edge_KDE_value += GetEdgeKDEValue(model, e);
        }

//    }   else if (model.mask_switch == Param::ON) {
//        for (int e : model.visit_edges) {
//            edge_KDE_value += GetEdgeKDEValue(model, e);
//        }
    }   else {
        for (int e = 0;e < model.m; e++) {
            edge_KDE_value += GetEdgeKDEValue(model, e);
        }
    }

    model.cur_lixel.KDE_value = edge_KDE_value;
}

void InitItemMasking(Model &model) {
    double length = model.edge_set[model.cur_lixel.edge_index]->length;
    double interval = model.lixel_interval_length;
    model.lixels_aggregation.clear();
    model.lixels_aggregation.resize((int)((length + interval / 2) / interval), 0);
    LixelAggregate(model, model.full_masked_c_edge, model.lixels_aggregation, 1);
    LixelAggregate(model, model.full_masked_d_edge, model.lixels_aggregation, 0);
    for (int i = 2;i < model.lixels_aggregation.size(); i++) {
        model.lixels_aggregation[i] += model.lixels_aggregation[i-1];
    }
    for (int i = 1;i < model.lixels_aggregation.size(); i++) {
        model.lixels_aggregation[i] += model.lixels_aggregation[i-1];
    }
//    for (int i = 0;i < model.lixels_aggregation.size(); i++) {
//        cout << "my answer " << model.lixels_aggregation[i] << endl;
//    }
}

double SimpleMaskKernelValue(Model &model, double interval, int &coe_c, int &coe_d) {
    return MaskKernelValue(model, coe_c * interval, coe_d * interval,
                           model.item1, model.item2, model.item3, model.item4);
}

void LixelAggregate(Model &model, vector<int> edges, vector<double> &vec, int is_node_c) {
    double length = model.edge_set[model.cur_lixel.edge_index]->length;
    double interval = model.lixel_interval_length;

    for (int i = 0;i < edges.size(); i++) {
        int e = edges[i];
        int node_index_c, coe_c, coe_d;
        if (is_node_c == 1) {
            node_index_c = model.edge_set[e]->node_1;
            coe_c = 1; coe_d = 0;
        }   else {
            node_index_c = model.edge_set[e]->node_2;
            coe_c = 0; coe_d = 1;
        }
        model.lixels_aggregation[0] += GetEdgeKDEValue(model, e);
        double dist_a_c = model.sp_node_a_vec[node_index_c].cur_sp_value;
        double dist_b_c = model.sp_node_b_vec[node_index_c].cur_sp_value;
        double dist_a_b = length;
        double breakpoint = (dist_a_c + dist_b_c + dist_a_b) / 2 - dist_a_c;
        if (breakpoint <= interval / 2) {
            vec[1] -= SimpleMaskKernelValue(model, interval, coe_c, coe_d);

        }   else if (breakpoint >= vec.size() * interval - interval / 2) {
            vec[1] += SimpleMaskKernelValue(model, interval, coe_c, coe_d);

        }   else {
            int in = (int)((breakpoint + interval / 2) /interval);
            double xx = breakpoint - (in * interval - interval / 2);
            double yy = interval - xx;
            vec[1] += SimpleMaskKernelValue(model, interval, coe_c, coe_d);
            if (in > 0) vec[in] += SimpleMaskKernelValue(model, xx - yy - interval, coe_c, coe_d);
            if (in + 1 < vec.size()) vec[in + 1] += SimpleMaskKernelValue(model, - interval - xx + yy, coe_c, coe_d);
        }
    }

}

void GetItemMasking(Model &model, int edge_index, Item &item1, Item &item2, Item &item3, Item &item4) {
    model.item1 = item1;
    model.item2 = item2;
    model.item3 = item3;
    model.item4 = item4;
}

double GetEdgeKDEValue(Model &model, int edge_index) {
    int node_index_c = model.edge_set[edge_index]->node_1;
    int node_index_d = model.edge_set[edge_index]->node_2;
    double dist_q_c = model.sp_lixel_vec[node_index_c].cur_sp_value;
    double dist_q_d = model.sp_lixel_vec[node_index_d].cur_sp_value;
    if (dist_q_c >= model.bandwidth && dist_q_d >= model.bandwidth) return 0;
    bool on_edge = (model.cur_lixel.edge_index == edge_index);

    if (model.method_type == Param::RQS || model.method_type == Param::SPS) {
        vector<pair<double, double>> q_p;
        model.edge_set[edge_index]->RQS_GetDistSum(model.start_time, model.end_time, dist_q_c, dist_q_d, q_p, on_edge);

        double edge_KA_value = 0;
        for (int i = 0;i < q_p.size(); i++) {
            if (q_p[i].first < model.bandwidth) {
                model.num++;
                edge_KA_value += KernelValue(model, q_p[i]);
            }
        }

        return edge_KA_value;
    }

    if (on_edge) {
        // dis1 <--bandwidth--> lixel=dis0 <--bandwidth--> dis2
        double dis0 = model.cur_lixel.dist_1;
        double dis1 = model.cur_lixel.dist_1 - model.bandwidth;
        double dis2 = model.cur_lixel.dist_1 + model.bandwidth;
        double middle_time = (model.start_time + model.end_time) / 2.0;
        Item item1 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, true, dis2, dis0 + EPS);
        Item item2 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, true, dis2, dis0 + EPS);
        Item item3 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, true, dis0, dis1);
        Item item4 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, true, dis0, dis1);
        return KernelValue(model, -dist_q_c, -dist_q_d, item1, item2, item3, item4);
//        return  KernelValue1(model, -dist_q_c, item1.parameters[2], item1.parameters[0], item1.parameters[4], item1.num) +
//                KernelValue2(model, -dist_q_c, item2.parameters[2], item2.parameters[0], item2.parameters[4], item2.num) +
//                KernelValue1(model, -dist_q_d, item3.parameters[3], item3.parameters[1], item3.parameters[4], item3.num) +
//                KernelValue2(model, -dist_q_d, item4.parameters[3], item4.parameters[1], item4.parameters[4], item4.num);
    }

    double breakpoint = (dist_q_c + model.edge_set[edge_index]->length + dist_q_d) / 2;
    double breakpoint_c = breakpoint - dist_q_c;
    double breakpoint_d = breakpoint - dist_q_d;
    double dist_c = min(breakpoint_c, model.bandwidth - dist_q_c);
    double dist_d = min(breakpoint_d, model.bandwidth - dist_q_d);

    double middle_time = (model.start_time + model.end_time) / 2.0;
//    Item item1 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, true, dist_c - EPS);
//    Item item2 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, true, dist_c - EPS);
//    Item item3 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, false, dist_d);
//    Item item4 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, false, dist_d);
    Item item1 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, true, dist_c);
    Item item2 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, true, dist_c);
    Item item3 = model.edge_set[edge_index]->GetDistSum(model.start_time, middle_time - EPS, false, dist_d - EPS);
    Item item4 = model.edge_set[edge_index]->GetDistSum(middle_time, model.end_time, false, dist_d - EPS);
    if (model.mask_switch == Param::ON && model.edge_state == Param::NEW) GetItemMasking(model, edge_index, item1, item2, item3, item4);
//    if (model.real_switch == Param::ON && model.edge_state == Param::NEW) GetItemMasking(model, edge_index, item1, item2, item3, item4);
    return KernelValue(model, dist_q_c, dist_q_d, item1, item2, item3, item4);
//    return  Triangular::KernelValue1(model, dist_q_c, item1.parameters[2], item1.parameters[0], item1.parameters[4], item1.num) +
//            Triangular::KernelValue2(model, dist_q_c, item2.parameters[2], item2.parameters[0], item2.parameters[4], item2.num) +
//            Triangular::KernelValue1(model, dist_q_d, item3.parameters[3], item3.parameters[1], item3.parameters[4], item3.num) +
//            Triangular::KernelValue2(model, dist_q_d, item4.parameters[3], item4.parameters[1], item4.parameters[4], item4.num);
}

double KernelValue(Model &model, pair<double, double> dist_time) {
    switch (model.kernel_type) {
        case Param::Triangular:
            return Triangular::KernelValue(model, dist_time);
            break;
        case Param::Exponential:
            return Exponential::KernelValue(model, dist_time);
            break;
        case Param::Cosine:
            return Cosine::KernelValue(model, dist_time);
            break;
        default:
            assert(false);
    }
}

double MaskKernelValue(Model &model, double dist_q_c, double dist_q_d, Item &item1, Item &item2, Item &item3, Item &item4) {
    switch (model.kernel_type) {
        case Param::Triangular:
            return  Triangular::MaskKernelValue1(model, dist_q_c, item1.parameters[2], item1.parameters[0], item1.parameters[4], item1.num) +
                    Triangular::MaskKernelValue2(model, dist_q_c, item2.parameters[2], item2.parameters[0], item2.parameters[4], item2.num) +
                    Triangular::MaskKernelValue1(model, dist_q_d, item3.parameters[3], item3.parameters[1], item3.parameters[4], item3.num) +
                    Triangular::MaskKernelValue2(model, dist_q_d, item4.parameters[3], item4.parameters[1], item4.parameters[4], item4.num);
            break;
//        case Param::Exponential:
//            return  Exponential::KernelValue1(model, dist_q_c, item1.parameters[0]) +
//                    Exponential::KernelValue2(model, dist_q_c, item2.parameters[1]) +
//                    Exponential::KernelValue1(model, dist_q_d, item3.parameters[2]) +
//                    Exponential::KernelValue2(model, dist_q_d, item4.parameters[3]);
//            break;
//        case Param::Cosine:
//            return  Cosine::KernelValue1(model, dist_q_c, item1.parameters[0], item1.parameters[1], item1.parameters[2], item1.parameters[3]) +
//                    Cosine::KernelValue2(model, dist_q_c, item2.parameters[0], item2.parameters[1], item2.parameters[2], item2.parameters[3]) +
//                    Cosine::KernelValue1(model, dist_q_d, item3.parameters[4], item3.parameters[5], item3.parameters[6], item3.parameters[7]) +
//                    Cosine::KernelValue2(model, dist_q_d, item4.parameters[4], item4.parameters[5], item4.parameters[6], item4.parameters[7]);
//            break;
        default:
            assert(false);
    }
}

double Triangular::MaskKernelValue1(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;
    double ans2 = ( q_c / b_s ) * sum_t / b_t;
    double ans4 = (q_c / b_s) * (1 - q_t / b_t) * num;
    return - ans2 - ans4;
}

double Triangular::MaskKernelValue2(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;

    double ans2 = ( q_c / b_s ) * sum_t / b_t;
    double ans4 = ( q_c / b_s ) * (1 + q_t / b_t) * num;
    return ans2 - ans4;
}

double KernelValue(Model &model, double dist_q_c, double dist_q_d, Item &item1, Item &item2, Item &item3, Item &item4) {
    switch (model.kernel_type) {
        case Param::Triangular:
            return  Triangular::KernelValue1(model, dist_q_c, item1.parameters[2], item1.parameters[0], item1.parameters[4], item1.num) +
                    Triangular::KernelValue2(model, dist_q_c, item2.parameters[2], item2.parameters[0], item2.parameters[4], item2.num) +
                    Triangular::KernelValue1(model, dist_q_d, item3.parameters[3], item3.parameters[1], item3.parameters[4], item3.num) +
                    Triangular::KernelValue2(model, dist_q_d, item4.parameters[3], item4.parameters[1], item4.parameters[4], item4.num);
            break;
        case Param::Exponential:
            return  Exponential::KernelValue1(model, dist_q_c, item1.parameters[0]) +
                    Exponential::KernelValue2(model, dist_q_c, item2.parameters[1]) +
                    Exponential::KernelValue1(model, dist_q_d, item3.parameters[2]) +
                    Exponential::KernelValue2(model, dist_q_d, item4.parameters[3]);
            break;
        case Param::Cosine:
            return  Cosine::KernelValue1(model, dist_q_c, item1.parameters[0], item1.parameters[1], item1.parameters[2], item1.parameters[3]) +
                    Cosine::KernelValue2(model, dist_q_c, item2.parameters[0], item2.parameters[1], item2.parameters[2], item2.parameters[3]) +
                    Cosine::KernelValue1(model, dist_q_d, item3.parameters[4], item3.parameters[5], item3.parameters[6], item3.parameters[7]) +
                    Cosine::KernelValue2(model, dist_q_d, item4.parameters[4], item4.parameters[5], item4.parameters[6], item4.parameters[7]);
            break;
        default:
            assert(false);
    }
}

double Triangular::KernelValue(Model &model, pair<double, double> dist_time) {
    double b_s = model.bandwidth;
    double b_t = (model.end_time - model.start_time) / 2.0;

    return (1 - dist_time.first / b_s) * (1 - dist_time.second / b_t);
}

double Triangular::KernelValue1(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;
    return (-1 * sum_ct + (b_s - q_c) * sum_t - (b_t - q_t) * sum_c + (b_s - q_c) * (b_t - q_t) * num) / b_s / b_t;
}

double Triangular::KernelValue2(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;

    double ans1 = sum_ct / b_s / b_t;
    double ans2 = (1 - q_c / b_s) * sum_t / b_t;
    double ans3 = (1 + q_t / b_t) * sum_c / b_s;
    double ans4 = (1 - q_c / b_s) * (1 + q_t / b_t) * num;
    return ans1 - ans2 - ans3 + ans4;
}

double Exponential::KernelValue(Model &model, pair<double, double> dist_time) {
    double b_s = model.bandwidth;
    double b_t = (model.end_time - model.start_time) / 2.0;
    return exp(-dist_time.first / b_s) * exp(-dist_time.second / b_t);
}

double Exponential::KernelValue1(Model &model, double q_c, double sum_exp) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;

    return sum_exp * exp(-q_c / b_s) * exp(-q_t / b_t);
}

double Exponential::KernelValue2(Model &model, double q_c, double sum_exp) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;

    return sum_exp * exp(-q_c / b_s) * exp(+q_t / b_t);
}

double Cosine::KernelValue(Model &model, pair<double, double> dist_time) {
    double b_s = model.bandwidth;
    double b_t = (model.end_time - model.start_time) / 2.0;
    return cos(dist_time.first / b_s) * cos(dist_time.second / b_t);
}

double Cosine::KernelValue1(Model &model, double q_c, double sum_cos_cos, double sum_sin_cos, double sum_cos_sin, double sum_sin_sin) {
    double b_s = model.bandwidth;
    double q_t = (model.start_time + model.end_time) / 2.0;
    double b_t = (model.end_time - model.start_time) / 2.0;

    double ans1 = cos(q_c / b_s) * cos(q_t / b_t) * sum_cos_cos;
    double ans2 = sin(q_c / b_s) * cos(q_t / b_t) * sum_sin_cos;
    double ans3 = cos(q_c / b_s) * sin(q_t / b_t) * sum_cos_sin;
    double ans4 = sin(q_c / b_s) * sin(q_t / b_t) * sum_sin_sin;
    return ans1 - ans2 + ans3 - ans4;
}

double Cosine::KernelValue2(Model &model, double q_c, double sum_cos_cos, double sum_sin_cos, double sum_cos_sin, double sum_sin_sin) {
    return KernelValue1(model, q_c, sum_cos_cos, sum_sin_cos, sum_cos_sin, sum_sin_sin);
}