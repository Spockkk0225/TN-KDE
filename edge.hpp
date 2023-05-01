#ifndef edge_hpp
#define edge_hpp

#include "library.hpp"

const int para_size = 8;

struct Item {
    int num = 0;
    double parameters[para_size] = {0};
//    vector<double> parameters{0,0,0,0,0};
    
    Item operator +(const Item &rhs) {
        Item ret;
        ret.num = this->num + rhs.num;
        for (int i = 0;i < para_size; i++) {
            ret.parameters[i] = this->parameters[i] + rhs.parameters[i];
        }
        return ret;
    }
    
    Item operator -(const Item &rhs) {
        Item ret;
        ret.num = this->num - rhs.num;
        for (int i = 0;i < para_size; i++) {
            ret.parameters[i] = this->parameters[i] - rhs.parameters[i];
        }
        return ret;
    }
    
    Item& operator +=(const Item &rhs) {
        this->num += rhs.num;
        for (int i = 0;i < para_size; i++) {
            this->parameters[i] += rhs.parameters[i];
        }
        return *this;
    }
    
    Item& operator -=(const Item &rhs) {
        this->num -= rhs.num;
        for (int i = 0;i < para_size; i++) {
            this->parameters[i] -= rhs.parameters[i];
        }
        return *this;
    }
};

struct Point {
    int edge_index;
    
    double dist_1, dist_2;
    double dist;
    double time;
    double b_s, b_t;
    
    Item item;
    
    void InitItem(Param::Kernel kernel_type) {
        item.num = 1;
        if (kernel_type == Param::Triangular) {
            item.parameters[0] = dist_1;
            item.parameters[1] = dist_2;
            item.parameters[2] = dist_1 * time;
            item.parameters[3] = dist_2 * time;
            item.parameters[4] = time;
            
        }   else if (kernel_type == Param::Exponential) {
            item.parameters[0] = exp(-dist_1 / b_s) * exp(+time / b_t);
            item.parameters[1] = exp(-dist_1 / b_s) * exp(-time / b_t);
            item.parameters[2] = exp(-dist_2 / b_s) * exp(+time / b_t);
            item.parameters[3] = exp(-dist_2 / b_s) * exp(-time / b_t);
            
        }   else if (kernel_type == Param::Cosine) {
            item.parameters[0] = cos(dist_1 / b_s) * cos(time / b_t);
            item.parameters[1] = sin(dist_1 / b_s) * cos(time / b_t);
            item.parameters[2] = cos(dist_1 / b_s) * sin(time / b_t);
            item.parameters[3] = sin(dist_1 / b_s) * sin(time / b_t);
            item.parameters[4] = cos(dist_2 / b_s) * cos(time / b_t);
            item.parameters[5] = sin(dist_2 / b_s) * cos(time / b_t);
            item.parameters[6] = cos(dist_2 / b_s) * sin(time / b_t);
            item.parameters[7] = sin(dist_2 / b_s) * sin(time / b_t);
            
        }   else {
            assert(false);
        }
    }
};
bool PointTimeCmp(const Point &x, const Point &y);
bool PointDistCmp(const Point &x, const Point &y);

struct Lixel : Point {
    double KDE_value;
};

struct SPNode {
    int node_index;
    double cur_sp_value;
};

// cmp for minimum heap of Dijkstra
struct SPNode_cmp {
    bool operator()(SPNode &n1, SPNode &n2) {
        return n1.cur_sp_value > n2.cur_sp_value;
    }
};

class Edge {
public:
    double test = 0;
    
    int node_1, node_2;
    double length;
    double b_s, b_t;
    vector<Point> point_set;
    
    int siz_t, siz_d;
    vector<double> time_axis, time_axis_inv;
    vector<double> dist_axis, dist_axis_inv;
    void InitAxis(Param::Kernel kernel_type);
    
    // Solution Baseline
    void RQS_GetDistSum(int start_time, int end_time, double dist_q_c, double dist_q_d, vector<pair<double, double>> &q_p, bool on_edge);

    double start_time, end_time;
    void SetTime(double the_start_time, double the_end_time);

    virtual void Init(Param::Kernel kernel_type) = 0;
    virtual Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1) = 0;
    
    // used for Dynamic Deep
    virtual void SetH(int h) {}

};

class Baseline : public Edge {
public:
    // useless, just implement the virtual function
    void Init(Param::Kernel kernel_type) { return; }
    Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1) { return Item(); }
};

class RAS : public Edge {
public:
    void Init(Param::Kernel kernel_type);
    Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1);
private:
    vector<vector<Item>> pre_sum;
    Item GetPreSum(int t1, int t2, int d1, int d2);
};

struct TreeDist {
    int l_child = -1, r_child = -1;
    double l_int = -1, r_int = -1;
    int tree_time_root = -1;
    vector<Point> time_ordered_point;
};
struct TreeTime {
    int l_child = -1, r_child = -1;
    double l_int = -1, r_int = -1;
    Item item;
};

class RTS : public Edge {
public:
    void Init(Param::Kernel kernel_type);
    Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1);
private:
    int tree_dist_root = -1;
    vector<TreeDist> tree_dist;
    vector<TreeTime> tree_time;
    Item GetSum(double ls, double rs, double lt, double rt);
    
    void MergeSort(vector<Point> &t, vector<Point> &s1, vector<Point> &s2);
    int SpatialConstruct(int index, int l, int r, vector<Point> &point_set);
    int TemporalConstruct(int index, int l, int r, vector<Point> &point_set);
    Item SpatialQuery(int index, double ls, double rs, double lt, double rt);
    Item TemporalQuery(int index, double lt, double rt);
};

class RFS : public Edge {
public:
    void Init(Param::Kernel kernel_type);
    Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1);
private:
    int timeline_index_L_1, timeline_index_R_1;
    int timeline_index_L_2, timeline_index_R_2;
    vector<int> timeline;
    vector<int> tree_root;
    vector<TreeTime> tree;
    Item GetSum(int L1, int R1, double L2, double R2);
    int Build(int l, int r);
    int Insert(int l, int r, int root_l, int L, Item &item);
    Item Query(int l, int r, int root_l, int root_r, double L, double R);
};

struct LazyDeep {
    double l, r;
    int root_l, root_r;
};

class DRFS : public Edge {
public:
    void Init(Param::Kernel kernel_type);
    Item GetDistSum(int start_time, int end_time, bool is_left, double dist_r, double dist_l = -1);
private:
    int H = 1;
    int oldH = 1;
    int Hmax = 1;
    int Hinsert = 1;
    double Dmin, Dmax;
    void SetH(int h);
    
    int timeline_index_L_1, timeline_index_R_1;
    int timeline_index_L_2, timeline_index_R_2;
    vector<int> timeline;
    vector<int> tree_root;
    vector<TreeTime> tree;
    vector<LazyDeep> lazy_deep;
    Item GetSum(int L1, int R1, double L2, double R2);
    int Build(int h, double l, double r);
    int Insert(int h, double l, double r, int root_l, double L, Item &item);
    Item Query(int h, double l, double r, int root_l, int root_r, double L, double R);
};

#endif /* edge_hpp */
