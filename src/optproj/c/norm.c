#include <math.h>
#include "norm.h"

double l2_norm(const double *arr, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double v = arr[i];
        s += v * v;
    }
    return sqrt(s);
}
