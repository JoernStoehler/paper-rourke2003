#include <math.h>
#include "numerics.h"

double l2_norm(const double *arr, int n) {
    double s = 0.0;
    for (int i = 0; i < n; ++i) {
        double v = arr[i];
        s += v * v;
    }
    return sqrt(s);
}

void inradius(
    // Polygon
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    // Outputs
    double *r, // inradius
    double drdx[], // gradient w.r.t. x, [s]
    double drdy[]  // gradient w.r.t. y, [s]
){
    /*
     * Algorithm for finding the inradius (and the incenter) of a convex polygon:
     * 
     */
}

void circumradius(
    // Polygon
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    // Outputs
    double *R, // circumradius
    double dRdx[], // gradient w.r.t. x, [s]
    double dRdy[]  // gradient w.r.t. y, [s]
) {}

void optimal_partition(
    // Vertices
    const int n,
    const double x[], // [n]
    const double y[], // [n]
    // Search Cutoffs
    const int max_m, // max number of polygons
    const int max_s, // max number of vertices in a polygon
    // Output
    int *m, // number of polygons
    double s[], // number of vertices in each polygon, [max_m]
    int p[], // indices of vertices in each polygon, [max_m*max_s] = [max_m][max_s]
    double *L, // total aspect ratio
    int *k // index of polygon with worst aspect ratio
) {}

void gradient_step(
    // Vertices
    const int n,
    double x[], // [n]
    double y[], // [n]
    // Search Cutoffs
    const int max_m, // max number of polygons
    const int max_s, // max number of vertices in a polygon
    // Gradient Descend Step
    const double alpha // step size
) {
    // 1. Find optimal partition, and its relevant polygon
    // 2. Compute gradient of aspect ratio of that polygon
    // 3. Take a step in the negative gradient direction

    int m; // number of polygons
    double s[max_m]; // number of vertices in each polygon, [max_m]
    int p[max_m * max_s]; // indices of vertices in each polygon, [max_m*max_s] = [max_m][max_s]
    double L; // total aspect ratio
    int k; // index of polygon with worst aspect ratio

    optimal_partition(n, x, y, max_m, max_s, &m, s, p, &L, &k);

    int sk = (int)s[k]; // number of vertices in polygon k
    int *pk = &p[k * max_s]; // indices of vertices in polygon k
    double xk[max_s]; // x coordinates of vertices in polygon k
    double yk[max_s]; // y coordinates of vertices in polygon k
    for (int i = 0; i < sk; i++) {
        xk[i] = x[pk[i]];
        yk[i] = y[pk[i]];
    }

    double r, R; // inradius and circumradius
    double drdx[max_s], drdy[max_s]; // gradient of inradius w.r.t. xk, yk
    double dRdx[max_s], dRdy[max_s]; // gradient of circumradius w.r.t. xk, yk
    inradius(sk, xk, yk, &r, drdx, drdy);
    circumradius(sk, xk, yk, &R, dRdx, dRdy);

    // A = R / r > 1
    // dA/dx = (r * dR/dx - R * dr/dx) / r^2
    // dA/dy = (r * dR/dy - R * dr/dy) / r^2
    double dAdx[max_s];
    double dAdy[max_s];
    for (int i = 0; i < sk; i++) {
        dAdx[i] = (r * dRdx[i] - R * drdx[i]) / (r * r);
        dAdy[i] = (r * dRdy[i] - R * drdy[i]) / (r * r);
    }

    // Take a step in the negative gradient direction
    for (int i = 0; i < sk; i++) {
        x[pk[i]] -= alpha * dAdx[i];
        y[pk[i]] -= alpha * dAdy[i];
    }

    return;
}