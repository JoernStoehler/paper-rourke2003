#ifndef OPTPROJ_NUMERICS_H
#define OPTPROJ_NUMERICS_H

double l2_norm(const double *arr, int n);

void inradius(
    // Polygon
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    // Outputs
    double *r, // inradius
    double drdx[], // gradient w.r.t. x, [s]
    double drdy[]  // gradient w.r.t. y, [s]
);

void circumradius(
    // Polygon
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    // Outputs
    double *R, // circumradius
    double dRdx[], // gradient w.r.t. x, [s]
    double dRdy[]  // gradient w.r.t. y, [s]
);

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
);

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
);

#endif