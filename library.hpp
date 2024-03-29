#ifndef library_hpp
#define library_hpp

#include <math.h>
#include <queue>
#include <string.h>
#include <assert.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

const int INF = 0x3f3f3f3f;
const double EPS = 1e-5;

namespace Param {
    enum Method{RQS, SPS, RAS, RTS, RFS, DRFS, ADA};
    enum Kernel{Triangular, Exponential, Cosine, Gaussian, Epanechnikov, Quartic};
    enum Prune{OFF, ON};
    enum State{OLD, NEW};
}

using namespace std;

#endif /* library_hpp */