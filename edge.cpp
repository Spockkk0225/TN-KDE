#include "edge.hpp"

void Edge::RQS_GetDistSum(int start_time, int end_time, double dist_q_c, double dist_q_d, vector<pair<double, double>> &q_p, bool on_edge) {
    q_p.clear();
    
    for (int i = 0; i < point_set.size(); i++) {
        if (point_set[i].time < start_time || point_set[i].time > end_time) continue;
        
        double dist;
        
        // special case: the lixel is on the edge
        if (on_edge == true) {
            dist = fabs(point_set[i].dist_1 - dist_q_c);
            
        }   else {
            dist = min(dist_q_c + point_set[i].dist_1,
                       dist_q_d + point_set[i].dist_2);
        }
        double time = fabs(point_set[i].time - (start_time + end_time) / 2.0);
        
        q_p.push_back(make_pair(dist, time));
    }
}

void Edge::SetTime(double the_start_time, double the_end_time) {
    this->start_time = the_start_time;
    this->end_time = the_end_time;
}

bool PointTimeCmp(const Point &x, const Point &y) { return x.time < y.time; }
bool PointDistCmp(const Point &x, const Point &y) { return x.dist < y.dist; }

void Edge::InitAxis(Param::Kernel kernel_type) {
    time_axis.clear();
    dist_axis.clear();
    for (int p = 0;p < point_set.size(); p++) {
        point_set[p].dist = point_set[p].dist_1;
        point_set[p].b_s = b_s;
        point_set[p].b_t = b_t;
        point_set[p].InitItem(kernel_type);
        time_axis.push_back(point_set[p].time);
        dist_axis.push_back(point_set[p].dist);
        dist_axis_inv.push_back(length - point_set[p].dist);
    }
    sort(time_axis.begin(), time_axis.end());
    sort(dist_axis.begin(), dist_axis.end());
    sort(dist_axis_inv.begin(), dist_axis_inv.end());
    vector<double>::iterator pos1 = unique(time_axis.begin(), time_axis.end());
    vector<double>::iterator pos2 = unique(dist_axis.begin(), dist_axis.end());
    vector<double>::iterator pos3 = unique(dist_axis_inv.begin(), dist_axis_inv.end());
    time_axis.erase(pos1, time_axis.end());
    dist_axis.erase(pos2, dist_axis.end());
    dist_axis_inv.erase(pos3, dist_axis_inv.end());
}

/* ----------------------------------------------------------- */

void RAS::Init(Param::Kernel kernel_type) {
    InitAxis(kernel_type);
    
    Item temp_item;
    vector<Item> temp_vec;
    temp_vec.resize(dist_axis.size(), temp_item);
    pre_sum.resize(time_axis.size(), temp_vec);
    for (int p = 0;p < point_set.size(); p++) {
        int time_index = (int)(lower_bound(time_axis.begin(), time_axis.end(), point_set[p].time) - time_axis.begin());
        int dist_index = (int)(lower_bound(dist_axis.begin(), dist_axis.end(), point_set[p].dist) - dist_axis.begin());
        pre_sum[time_index][dist_index] += point_set[p].item;
    }
    
    for (int i = 0;i < time_axis.size(); i++) {
        for (int j = 0;j < dist_axis.size(); j++) {
            if (i > 0) pre_sum[i][j] += pre_sum[i-1][j];
            if (j > 0) pre_sum[i][j] += pre_sum[i][j-1];
            if (i > 0 && j > 0) pre_sum[i][j] -= pre_sum[i-1][j-1];
        }
    }
}

Item RAS::GetPreSum(int t1, int t2, int d1, int d2) {
    Item ret;
    if (t2 < t1 || d2 < d1) return ret;
    // pre_sum[t2][d2] - pre_sum[t2][d1-1] - pre_sum[t1-1][d2] + pre_sum[t1-1][d1-1]
    ret = pre_sum[t2][d2];
    if (d1 > 0) ret -= pre_sum[t2][d1-1];
    if (t1 > 0) ret -= pre_sum[t1-1][d2];
    if (d1 > 0 && t1 > 0) ret += pre_sum[t1-1][d1-1];
    return ret;
}

Item RAS::GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l) {
    Item ret;
    // if (time_axis.size() == 0) return Item();
    int time_index_1 = (int)(lower_bound(time_axis.begin(), time_axis.end(), start_time) - time_axis.begin());
    int time_index_2 = (int)(upper_bound(time_axis.begin(), time_axis.end(), end_time) - time_axis.begin()) - 1;

    if (is_left) {
        int dist_index_1 = (int)(lower_bound(dist_axis.begin(), dist_axis.end(), dist_l) - dist_axis.begin());
        int dist_index_2 = (int)(upper_bound(dist_axis.begin(), dist_axis.end(), dist_r) - dist_axis.begin()) - 1;
        ret = GetPreSum(time_index_1, time_index_2, dist_index_1, dist_index_2);
        return ret;
        
    }   else {
        
        int dist_index_1 = (int)(lower_bound(dist_axis_inv.begin(), dist_axis_inv.end(), dist_l) - dist_axis_inv.begin());
        int dist_index_2 = (int)(upper_bound(dist_axis_inv.begin(), dist_axis_inv.end(), dist_r) - dist_axis_inv.begin()) - 1;
        dist_index_1 = (int)dist_axis_inv.size() - dist_index_1 - 1;
        dist_index_2 = (int)dist_axis_inv.size() - dist_index_2 - 1;
        ret = GetPreSum(time_index_1, time_index_2, dist_index_2, dist_index_1);
        return ret;
    }
}

/* ----------------------------------------------------------- */

void RTS::Init(Param::Kernel kernel_type) {
    InitAxis(kernel_type);
    
    test = test + point_set.size();
    sort(point_set.begin(), point_set.end(), PointDistCmp);
    tree_dist_root = SpatialConstruct(tree_dist_root, 0, (int)point_set.size() - 1, point_set);
    
//    for (int p = 0;p < point_set.size(); p++) {
//        int time_index = (int)(lower_bound(time_axis.begin(), time_axis.end(), point_set[p].time) - time_axis.begin());
//        int dist_index = (int)(lower_bound(dist_axis.begin(), dist_axis.end(), point_set[p].dist) - dist_axis.begin());
//        tree_time_root = InsertTreeTime(tree_time_root, 0, (int)time_axis.size()-1, point_set[p].item, time_index, dist_index);
//    }
}

Item RTS::GetSum(double ls, double rs, double lt, double rt) {
    if (ls > rs || lt > rt) return Item();
    return SpatialQuery(tree_dist_root, ls, rs, lt, rt);
}

Item RTS::GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l) {
    if (is_left) {
        return GetSum(dist_l, dist_r, start_time, end_time);
    }   else {
        return GetSum(length - dist_r, length - dist_l, start_time, end_time);
    }
}

void RTS::MergeSort(vector<Point> &t, vector<Point> &s1, vector<Point> &s2) {
    t.clear();
    int p1 = 0, p2 = 0;
    while (p1 < s1.size() && p2 < s2.size()) {
        if (s1[p1].time < s2[p2].time) {
            t.push_back(s1[p1]);
            p1++;
        }   else {
            t.push_back(s2[p2]);
            p2++;
        }
    }
    while (p1 < s1.size()) {
        t.push_back(s1[p1]);
        p1++;
    }
    while (p2 < s2.size()) {
        t.push_back(s2[p2]);
        p2++;
    }
}

// do not use &(reference)! vector.push_back will change the address
// use return value to save the index
int RTS::SpatialConstruct(int index, int l, int r, vector<Point> &point_set) {
    // 7 2 3
    if (l > r) return -1;
    if (index == -1) {
        index = (int) tree_dist.size();
        TreeDist temp_tree_dist;
        tree_dist.push_back(temp_tree_dist);
    }
    tree_dist[index].l_int = point_set[l].dist;
    tree_dist[index].r_int = point_set[r].dist;
  //  cout << index << " " << l << " " << r << endl;
    
    if (l != r) {
        int mid = (l + r) >> 1;
        int lc = SpatialConstruct(tree_dist[index].l_child, l, mid, point_set);
        tree_dist[index].l_child = lc;
        int rc = SpatialConstruct(tree_dist[index].r_child, mid + 1, r, point_set);
        tree_dist[index].r_child = rc;
        //cout << "**index = " << index << " " << l << " " << r << endl;
        //cout << tree_dist[tree_dist[index].l_child].time_ordered_point.size() << endl;
        //cout << tree_dist[tree_dist[index].r_child].time_ordered_point.size() << endl;
        MergeSort(tree_dist[index].time_ordered_point, tree_dist[tree_dist[index].l_child].time_ordered_point, tree_dist[tree_dist[index].r_child].time_ordered_point);
        //cout << "--index = " << index << endl;
    }   else {
        tree_dist[index].time_ordered_point.clear();
        tree_dist[index].time_ordered_point.push_back(point_set[l]);
    }
    int rt = TemporalConstruct(tree_dist[index].tree_time_root, 0, (int)tree_dist[index].time_ordered_point.size()-1, tree_dist[index].time_ordered_point);
    tree_dist[index].tree_time_root = rt;
    return index;
}

int RTS::TemporalConstruct(int index, int l, int r, vector<Point> &point_set) {
    if (l > r) return -1;
    if (index == -1) {
        index = (int)tree_time.size();
        TreeTime temp_tree_time;
        tree_time.push_back(temp_tree_time);
    }
    tree_time[index].l_int = point_set[l].time;
    tree_time[index].r_int = point_set[r].time;
    
    if (l == r) {
        tree_time[index].item = point_set[l].item;
    }   else {
        int mid = (l + r) >> 1;
        int lc = TemporalConstruct(tree_time[index].l_child, l, mid, point_set);
        tree_time[index].l_child = lc;
        int rc = TemporalConstruct(tree_time[index].r_child, mid + 1, r, point_set);
        tree_time[index].r_child = rc;
        tree_time[index].item = tree_time[tree_time[index].l_child].item + tree_time[tree_time[index].r_child].item;
    }
    return index;
}

Item RTS::SpatialQuery(int index, double ls, double rs, double lt, double rt) {
    if (ls > rs || lt > rt) return Item();
    if (index == -1) return Item();
    if (tree_dist[index].l_int > rs || tree_dist[index].r_int < ls) return Item();
    if (ls <= tree_dist[index].l_int && tree_dist[index].r_int <= rs) {
        return TemporalQuery(tree_dist[index].tree_time_root, lt, rt);
    }
    
    Item ret;
    if (ls <= tree_dist[tree_dist[index].l_child].r_int) ret += SpatialQuery(tree_dist[index].l_child, ls, rs, lt, rt);
    if (rs >= tree_dist[tree_dist[index].r_child].l_int) ret += SpatialQuery(tree_dist[index].r_child, ls, rs, lt, rt);
    return ret;
}

Item RTS::TemporalQuery(int index, double lt, double rt) {
    if (lt > rt) return Item();
    if (index == -1) return Item();
    if (tree_time[index].l_int > rt || tree_time[index].r_int < lt) return Item();
    if (lt <= tree_time[index].l_int && tree_time[index].r_int <= rt) {
        return tree_time[index].item;
    }
    
    Item ret;
    if (lt <= tree_time[tree_time[index].l_child].r_int) ret += TemporalQuery(tree_time[index].l_child, lt, rt);
    if (rt >= tree_time[tree_time[index].r_child].l_int) ret += TemporalQuery(tree_time[index].r_child, lt, rt);
    return ret;
}

/* ----------------------------------------------------------- */

void RFS::Init(Param::Kernel kernel_type) {
    sort(point_set.begin(), point_set.end(), PointTimeCmp);
    if (point_set.size() == 0) return;
    InitAxis(kernel_type);
    
    // shift one index
    // tree_root[0] is the root of a blank tree
    // point i's root is tree_root[i+1]
    tree_root.resize(point_set.size() + 1);
    timeline.resize(point_set.size() + 1);
    tree_root[0] = Build(0, (int)dist_axis.size()-1);
    timeline[0] = -1;
    for (int p = 0;p < point_set.size(); p++) {
        int dist_index = (int)(lower_bound(dist_axis.begin(), dist_axis.end(), point_set[p].dist) - dist_axis.begin());
        tree_root[p+1] = Insert(0, (int)dist_axis.size()-1, tree_root[p], dist_index, point_set[p].item);
        timeline[p+1] = point_set[p].time;
    }

    double middle_time = (start_time + end_time) / 2;
    timeline_index_L_1 = (int)(lower_bound(timeline.begin(), timeline.end(), start_time) - timeline.begin());
    timeline_index_R_1 = (int)(upper_bound(timeline.begin(), timeline.end(), middle_time - EPS) - timeline.begin()) - 1;
    timeline_index_L_2 = (int)(lower_bound(timeline.begin(), timeline.end(), middle_time) - timeline.begin());
    timeline_index_R_2 = (int)(upper_bound(timeline.begin(), timeline.end(), end_time) - timeline.begin()) - 1;
}

Item RFS::GetSum(int L1, int R1, double L2, double R2) {
    if (L1 > R1 || L2 > R2) return Item();
    if (R2 < dist_axis[0]) return Item();
    if (L2 > dist_axis[dist_axis.size()-1]) return Item();
    return Query(0, (int)dist_axis.size()-1, tree_root[L1-1], tree_root[R1], L2, R2);
}

Item RFS::GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l) {
    Item ret;
    if (timeline.size() == 0) return Item();
    
    int time_index_1, time_index_2;
    if (start_time == this->start_time) {
        time_index_1 = timeline_index_L_1;
        time_index_2 = timeline_index_R_1;
    }   else {
        time_index_1 = timeline_index_L_2;
        time_index_2 = timeline_index_R_2;
    }

    if (is_left) {
    //    int dist_index_1 = (int)(lower_bound(dist_axis.begin(), dist_axis.end(), dist_l) - dist_axis.begin());
    //    int dist_index_2 = (int)(upper_bound(dist_axis.begin(), dist_axis.end(), dist_r) - dist_axis.begin()) - 1;
        ret = GetSum(time_index_1, time_index_2, dist_l, dist_r);
        return ret;
        
    }   else {
        
    //    int dist_index_1 = (int)(lower_bound(dist_axis_inv.begin(), dist_axis_inv.end(), dist_l) - dist_axis_inv.begin());
    //    int dist_index_2 = (int)(upper_bound(dist_axis_inv.begin(), dist_axis_inv.end(), dist_r) - dist_axis_inv.begin()) - 1;
    //    dist_index_1 = (int)dist_axis_inv.size() - dist_index_1 - 1;
    //    dist_index_2 = (int)dist_axis_inv.size() - dist_index_2 - 1;
        ret = GetSum(time_index_1, time_index_2, length - dist_r, length - dist_l);
        return ret;
    }
}

int RFS::Build(int l, int r) {
    int root = (int)tree.size();
    TreeTime temp_tree_time;
    tree.push_back(temp_tree_time);
    return tree[root].l_child = tree[root].r_child = root;
//    if (l == r) return root;
//    int mid = (l + r) >> 1;
//    int lc = Build(l, mid);
//    tree[root].l_child = lc;
//    int rc = Build(mid+1, r);
//    tree[root].r_child = rc;
//    return root;
}

int RFS::Insert(int l, int r, int root_l, int L, Item &item) {
    int root_r = (int)tree.size();
    TreeTime temp_tree_time;
    tree.push_back(temp_tree_time);
    tree[root_r].l_child = tree[root_l].l_child;
    tree[root_r].r_child = tree[root_l].r_child;
    tree[root_r].item = tree[root_l].item + item;
    if (l == r) return root_r;
    int mid = (l + r) >> 1;
    if (L <= mid) {
        int lc = Insert(l, mid, tree[root_l].l_child, L, item);
        tree[root_r].l_child = lc;
    }   else {
        int rc = Insert(mid+1, r, tree[root_l].r_child, L, item);
        tree[root_r].r_child = rc;
    }
    return root_r;
}

Item RFS::Query(int l, int r, int root_l, int root_r, double L, double R) {
    if (L <= dist_axis[l] && dist_axis[r] <= R) return tree[root_r].item - tree[root_l].item;
    int mid = (l + r) >> 1;
    Item ret;
    if (L <= dist_axis[mid]) {
        ret += Query(l, mid, tree[root_l].l_child, tree[root_r].l_child, L, R);
    }
    if (dist_axis[mid+1] <= R) {
        ret += Query(mid+1, r, tree[root_l].r_child, tree[root_r].r_child, L, R);
    }
    return ret;
}

/* ------------------------------------- */

void DRFS::Init(Param::Kernel kernel_type) {
    sort(point_set.begin(), point_set.end(), PointTimeCmp);
    if (point_set.size() == 0) return;
    InitAxis(kernel_type);
    
    // shift one index
    // tree_root[0] is the root of a blank tree
    // point i's root is tree_root[i+1]
    // Dmin and Dmax is the full interval
    tree_root.resize(point_set.size() + 1);
    timeline.resize(point_set.size() + 1);
    Dmin = dist_axis[0];
    Dmax = dist_axis[(int)dist_axis.size()-1];
    tree_root[0] = Build(0, Dmin, Dmax);
    timeline[0] = -1;
    for (int p = 0;p < point_set.size(); p++) {
        tree_root[p+1] = Insert(0, Dmin, Dmax, tree_root[p], point_set[p].dist, point_set[p].item);
        timeline[p+1] = point_set[p].time;
    }
    Hinsert = 0;
    
    double middle_time = (start_time + end_time) / 2;
    timeline_index_L_1 = (int)(lower_bound(timeline.begin(), timeline.end(), start_time) - timeline.begin());
    timeline_index_R_1 = (int)(upper_bound(timeline.begin(), timeline.end(), middle_time - EPS) - timeline.begin()) - 1;
    timeline_index_L_2 = (int)(lower_bound(timeline.begin(), timeline.end(), middle_time) - timeline.begin());
    timeline_index_R_2 = (int)(upper_bound(timeline.begin(), timeline.end(), end_time) - timeline.begin()) - 1;
}

Item DRFS::GetSum(int L1, int R1, double L2, double R2) {
    if (L1 > R1 || L2 > R2) return Item();
    if (R2 < Dmin) return Item();
    if (L2 > Dmax) return Item();
    return Query(0, Dmin, Dmax, tree_root[L1-1], tree_root[R1], L2, R2);
}

Item DRFS::GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l) {
    Item ret;
    if (timeline.size() == 0) return Item();
    
    int time_index_1, time_index_2;
    if (start_time == this->start_time) {
        time_index_1 = timeline_index_L_1;
        time_index_2 = timeline_index_R_1;
    }   else {
        time_index_1 = timeline_index_L_2;
        time_index_2 = timeline_index_R_2;
    }

    if (is_left) {
        ret = GetSum(time_index_1, time_index_2, dist_l, dist_r);
        return ret;
        
    }   else {
        ret = GetSum(time_index_1, time_index_2, length - dist_r, length - dist_l);
        return ret;
    }
}

int DRFS::Build(int h, double l, double r) {
    int root = (int)tree.size();
    TreeTime temp_tree_time;
    tree.push_back(temp_tree_time);
    tree[root].l_child = tree[root].r_child = root;
    return root;
//    if (h > H) return root;
//    double mid = (l + r) / 2;
//    tree[root].l = Build(h+1, l, mid);
//    tree[root].r = Build(h+1, mid, r);
//    return root;
}

int DRFS::Insert(int h, double l, double r, int root_l, double L, Item &item) {
    int root_r = (int)tree.size();
    TreeTime temp_tree_time;
    tree.push_back(temp_tree_time);
    tree[root_r].l_child = tree[root_l].l_child;
    tree[root_r].r_child = tree[root_l].r_child;
    tree[root_r].item = tree[root_l].item + item;
    if (h == H) {
        if (Hinsert == 1) {
            lazy_deep.push_back(LazyDeep{l, r, root_l, root_r});
        }
        return root_r;
    }
    double mid = (l + r) / 2;
    if (L <= mid) {
        int lc = Insert(h+1, l, mid, tree[root_l].l_child, L, item);
        tree[root_r].l_child = lc;
    }   else {
        int rc = Insert(h+1, mid, r, tree[root_l].r_child, L, item);
        tree[root_r].r_child = rc;
    }
    return root_r;
}

Item DRFS::Query(int h, double l, double r, int root_l, int root_r, double L, double R) {
    if (L <= l && r <= R) return tree[root_r].item - tree[root_l].item;
    if (h == H) return Item();
    double mid = (l + r) / 2;
    Item ret;
    if (L <= mid) {
        ret += Query(h+1, l, mid, tree[root_l].l_child, tree[root_r].l_child, L, R);
    }
    if (mid < R) {
        ret += Query(h+1, mid, r, tree[root_l].r_child, tree[root_r].r_child, L, R);
    }
    return ret;
}

void DRFS::SetH(int h) {
    if (h <= Hmax) {
        H = h;
    }   else {
        Hmax = H = h;
        for (int p = 0;p < point_set.size(); p++) {
            double l = lazy_deep[p].l, r = lazy_deep[p].r;
            double mid = (l + r) / 2;
            double L = point_set[p].dist;
            int root_l = lazy_deep[p].root_l;
            int root_r = lazy_deep[p].root_r;
            tree[root_r].l_child = tree[root_l].l_child;
            tree[root_r].r_child = tree[root_l].r_child;
            if (L <= mid) {
                int lc = Insert(oldH+1, l, mid, tree[root_l].l_child, L, point_set[p].item);
                tree[root_r].l_child = lc;
            }   else {
                int rc = Insert(oldH+1, mid, r, tree[root_l].r_child, L, point_set[p].item);
                tree[root_r].r_child = rc;
            }
        }
        oldH = h;
    }
}
