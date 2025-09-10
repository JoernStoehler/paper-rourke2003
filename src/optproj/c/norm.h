#ifndef OPTPROJ_NORM_H
#define OPTPROJ_NORM_H

double l2_norm(const double *arr, int n);

void aspect_ratio_convex_polygon(
    // Vertices
    const int n,
    const double x[], // [n]
    const double y[], // [n]

    // Polygon
    const int s, // number of vertices in polygon
    const int p[], // indices of vertices, [n]

    double *L, // aspect ratio, output
    double dLdx[], // gradient w.r.t. x, output, [s]
    double dLdy[]  // gradient w.r.t. y, output, [s]
);

void aspect_ratio_partition(
    const int n,
    const double x[], // [n]
    const double y[], // [n]

    const int m, // number of polygons
    const int max_s, // max number of vertices in a polygon
    const int p[], // indices of vertices, [m*max_s]

)

#endif
