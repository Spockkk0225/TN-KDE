#ifndef kernel_hpp
#define kernel_hpp

#include "init.hpp"

void GetKDE(Model &model);
double GetEdgeKDEValue(Model &model, int edge_index);

double KernelValue(Model &model, pair<double, double> dist_time);
double KernelValue(Model &model, double dist_q_c, double dist_q_d, Item &item1, Item &item2, Item &item3, Item &item4);


namespace Triangular {
    double KernelValue(Model &model, pair<double, double> dist_time);
    double KernelValue1(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num);
    double KernelValue2(Model &model, double q_c, double sum_ct, double sum_c, double sum_t, int num);
};

namespace Exponential {
    double KernelValue(Model &model, pair<double, double> dist_time);
    double KernelValue1(Model &model, double q_c, double sum_exp);
    double KernelValue2(Model &model, double q_c, double sum_exp);
};

namespace Cosine {
    double KernelValue(Model &model, pair<double, double> dist_time);
    double KernelValue1(Model &model, double q_c, double sum_cos_cos, double sum_sin_cos, double sum_cos_sin, double sum_sin_sin);
    double KernelValue2(Model &model, double q_c, double sum_cos_cos, double sum_sin_cos, double sum_cos_sin, double sum_sin_sin);
};
#endif /* kernel_hpp */
