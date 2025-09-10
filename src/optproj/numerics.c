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

// CIRCUMRADIUS
// no derivatives

void circumradius_triangle(
    // inputs
    const double x[3],
    const double y[3],
    // outputs
    double *cx,
    double *cy,
    double *R
) {
    double d[3];
    d[0] = sqrt(pow(x[1] - x[2], 2) + pow(y[1] - y[2], 2));
    d[1] = sqrt(pow(x[0] - x[2], 2) + pow(y[0] - y[2], 2));
    d[2] = sqrt(pow(x[0] - x[1], 2) + pow(y[0] - y[1], 2));
    double dmax = fmax(fmax(d[0], d[1]), d[2]);
    double idmax = (dmax == d[0]) ? 0 : (dmax == d[1]) ? 1 : 2;

    // Case: obtuse triangle, use midpoint of longest edge d[i], which is opposite vertex i
    if (2*pow(d[idmax],2) > pow(d[0],2) + pow(d[1],2) + pow(d[2],2)) {
        int j = (idmax + 1) % 3;
        int k = (idmax + 2) % 3;
        *cx = 0.5 * (x[j] + x[k]);
        *cy = 0.5 * (y[j] + y[k]);
        *R = 0.5 * d[idmax];
        return;
    }

    // Case: acute or right triangle. points may be approximately degenerate so watch numerical stability
    // formula: https://en.wikipedia.org/wiki/Circumcircle
    
    // Cartesian coordinates for the circumcenter
    double D = 2.0 * (x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
    if (fabs(D) < 1e-12) {
        // points are approximately collinear, we approximate by the midpoint of the longest edge
        int j = (idmax + 1) % 3;
        int k = (idmax + 2) % 3;
        *cx = 0.5 * (x[j] + x[k]);
        *cy = 0.5 * (y[j] + y[k]);
        *R = 0.5 * d[idmax];
        return;
    }
    *cx = ( (pow(x[0], 2) + pow(y[0], 2)) * (y[1] - y[2]) 
           +(pow(x[1], 2) + pow(y[1], 2)) * (y[2] - y[0])
           +(pow(x[2], 2) + pow(y[2], 2)) * (y[0] - y[1])) / D;
    *cy = ( (pow(x[0], 2) + pow(y[0], 2)) * (x[2] - x[1]) 
           +(pow(x[1], 2) + pow(y[1], 2)) * (x[0] - x[2])
           +(pow(x[2], 2) + pow(y[2], 2)) * (x[1] - x[0])) / D;
    
    // Circumradius from pythagoras
    *R = sqrt(pow(x[0] - *cx, 2) + pow(y[0] - *cy, 2));
    return;
}

void circumradius_polygon(
    // inputs
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    const double atol // calculate circumradius only up to absolute tolerance atol
    // outputs
    double *cx, // circumcenter x
    double *cy, // circumcenter y
    double *R // circumradius
) {
    // We use Welzl's algorithm on CPU to find the smallest enclosing disk of the vertices, which will then enclose their convex hull as well
    // We randomize, since Welzl is only in expectation linear time, and worst case quadratic time
    // For simplicity we only randomize once at initialization, not at each recursion level
    int nP = s;
    double Px[s];
    double Py[s];
    int nR = 0;
    double Rx[0];
    double Ry[0];

    // shuffle
    {
        int indices[s];
        for (int i = 0; i < s; i++) {
            indices[i] = i;
        }
        for (int i = s - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            int temp = indices[i];
            indices[i] = indices[j];
            indices[j] = temp;
        }
        for (int i = 0; i < s; i++) {
            Px[i] = x[indices[i]];
            Py[i] = y[indices[i]];
        }
    }

    // Welzl's algorithm
    welzl(s, Px, Py, R, Rx, Ry, atol, cx, cy, R);

    // note: due to floating point math, and due to atol, we are no longer identifying the *unique* smallest enclosing disk
    // but just any of the disks that have similarly large radius as the exact smallest enclosing disk
    return;
}

void welzl(
    // inputs
    const int nP, // remaining points
    const double Px[], // [nP]
    const double Py[], // [nP]
    const int nR, // boundary points
    const double Rx[], // [nR]
    const double Ry[], // [nR]
    const double atol, // absolute tolerance for radius
    // outputs
    double *cx, // center x
    double *cy, // center y
    double *R // radius
){
    if (nR > 3) {
        error("welzl: nR > 3");
        return;
    } else if (nR == 3) { // trivial case, triangle
        circumradius_triangle(&Rx[0], &Ry[0], cx, cy, R);
        return;
    } else if (nP == 0) { // 0, 1, or 2 boundary points
        if (nR == 0) { // no points => zero radius disk
            *cx = 0.0;
            *cy = 0.0;
            *R = 0.0;
        } else if (nR == 1) { // one point => zero radius disk at that point
            *cx = Rx[0];
            *cy = Ry[0];
            *R = 0.0;
        } else if (nR == 2) { // two points => disk at midpoint with radius half distance
            *cx = 0.5 * (Rx[0] + Rx[1]);
            *cy = 0.5 * (Ry[0] + Ry[1]);
            *R = 0.5 * sqrt(pow(Rx[0] - Rx[1], 2) + pow(Ry[0] - Ry[1], 2));
        }
    } else {
        // pick a random point from P. We just use the last one for simplicity
        double px = Px[nP - 1];
        double py = Py[nP - 1];
        // recurse without p
        welzl(nP - 1, Px, Py, nR, Rx, Ry, atol, cx, cy, R);
        // if p is inside the circle, return
        if (pow(px - *cx, 2) + pow(py - *cy, 2) <= pow(*R + atol, 2)) {
            return;
        }
        // otherwise, p is on the boundary
        double Rx2[nR + 1];
        double Ry2[nR + 1];
        memcpy(Rx2, Rx, nR * sizeof(double));
        memcpy(Ry2, Ry, nR * sizeof(double));
        Rx2[nR] = px;
        Ry2[nR] = py;
        welzl(nP - 1, Px, Py, nR + 1, Rx2, Ry2, atol, cx, cy, R);
        return;
    }
}

// INRADIUS
// no derivatives

void inradius_polygon(
    // inputs
    const int s, // number of vertices in polygon
    const double x[], // [s]
    const double y[], // [s]
    const double atol, // calculate inradius only up to absolute tolerance atol
    // outputs
    double *cx, // incenter x
    double *cy, // incenter y
    double *r, // inradius
) {
    // Facts
    // - the medial axis of a convex polygon is a tree of straight line segments, which are defined by 2 edges of the polygon that they are equidistant to
    // - the tree has (graph theory) n leaves and n-2 internal nodes of degree 3
    // - if we track the minimal distance to the polygon boundary along the medial axis, then this function is continuous and piecewise affine
    // - along line segments, the function is still affine, and thus monotonic
    // - furthermore, since the polygon is convex, the function has no local maxima, but only global maxima
    // - as a result, the function only has global maxima, which is the center of the maximal inscribed circle (up to equally high, connected maxima)
    // - so the distance function starts at 0 at the polygon vertices, then increases linearly along the medial axis branches, changing slope at the internal nodes, until eventually the global maximum is reached
    // 

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
    double *cx, // circumcenter x
    double *cy, // circumcenter y
    double dRdx[], // gradient w.r.t. x, [s]
    double dRdy[]  // gradient w.r.t. y, [s]
) {
    /*
     * Assume the polygon is convex, s >= 3
     * We use Welzl's algorithm on CPU to find the smallest radius enclosing disk.
     */
    
}

// Welzl's algorithm for smallest enclosing disk

inline int argmax3d(const double a, const double b, const double c) {
    if (a >= b && a >= c) return 0;
    if (b >= a && b >= c) return 1;
    return 2;
}

void circumradius_triangle_obtuse(
    // inputs
    const double v[6], // [2][3], vertices of the triangle
    const int index_obtuse, // index of the obtuse angle vertex
    // outputs
    double c[2], // [2], circumcenter
    double *R, // radius
    double dRdv[6], // [2][3], gradients: circumradius w.r.t. vertices
) {
    const double *x = &v[0]; // x-coordinates
    const double *y = &v[3]; // y-coordinates

    int i = index_obtuse;
    int j = (i + 1) % 3;
    int k = (i + 2) % 3;

    // obtuse triangle, use midpoint of longest edge d[i], which is opposite vertex i
    c[0] = 0.5 * (x[j] + x[k]);
    c[1] = 0.5 * (y[j] + y[k]);
    double d = sqrt(pow(x[j] - x[k], 2) + pow(y[j] - y[k], 2));
    *R = 0.5 * d;

    // derivatives
    dRdv[i] = 0.0;
    dRdv[i + 3] = 0.0;
    dRdv[j] = 0.5 * (x[j] - x[k]) / d;
    dRdv[j + 3] = 0.5 * (y[j] - y[k]) / d;
    dRdv[k] = 0.5 * (x[k] - x[j]) / d;
    dRdv[k + 3] = 0.5 * (y[k] - y[j]) / d;

    // memset(dcdv, 0, 12 * sizeof(double));
    // // dcdv[:,:,i] = 0
    // // dcdv[:,:,j] = identity[2,2] * 0.5
    // // dcdv[:,:,k] = identity[2,2] * 0.5
    // dcdv[0 * 6 + 0 * 3 + j] = 0.5;
    // dcdv[1 * 6 + 1 * 3 + j] = 0.5;
    // dcdv[0 * 6 + 0 * 3 + k] = 0.5;
    // dcdv[1 * 6 + 1 * 3 + k] = 0.5;
    
    return;
}

void circumradius_triangle(
    const double v[6], // [2][3], inputs: vertices of the triangle
    double c[2], // [2], outputs: circumcenter
    double *R, // output: circumradius
    double dRdv[6], // [2][3], gradients: circumradius w.r.t. vertices
    double dcdv[12] // [2][2][3], gradients: circumcenter w.r.t. vertices (not calculated yet)
) {
    // Wikipedia: https://en.wikipedia.org/wiki/Circumcircle
    // 1. We check if the triangle is obtuse, in which case the circumcircle is the midpoint of the longest edge
    // 2. Otherwise, we use the circumcircle formula in barycentric coordinates

    const double *x = &v[0]; // x-coordinates
    const double *y = &v[3]; // y-coordinates

    double d[3];
    d[0] = sqrt(pow(x[1] - x[2], 2) + pow(y[1] - y[2], 2));
    d[1] = sqrt(pow(x[0] - x[2], 2) + pow(y[0] - y[2], 2));
    d[2] = sqrt(pow(x[0] - x[1], 2) + pow(y[0] - y[1], 2));

    // Test for obtuse triangle
    int argd = argmax3d(d[0], d[1], d[2]);
    if (pow(d[argd], 2) > pow(d[(argd + 1) % 3], 2) + pow(d[(argd + 2) % 3], 2)) {
        // obtuse triangle, use midpoint of longest edge d[i], which is opposite vertex i
        circumradius_triangle_obtuse(v, argd, &c[0], R, &dRdv[0], &dcdv[0]);
        return;
    }

    // We use a simple formula for cartesian coordinates for the circumcenter and circumradius
    double D = 2.0 * (x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
    if (fabs(D) < 1e-12) {
        // points are approximately collinear, we approximate by the midpoint of the longest edge
        circumradius_triangle_obtuse(v, argd, &c[0], R, &dRdv[0], &dcdv[0]);
        return;
    }

    // Non-obtuse and non-degenerate triangle, use the cartesian formula
    c[0] = ((pow(x[0], 2) + pow(y[0], 2)) * (y[1] - y[2]) 
           +(pow(x[1], 2) + pow(y[1], 2)) * (y[2] - y[0])
           +(pow(x[2], 2) + pow(y[2], 2)) * (y[0] - y[1])) / D;
    c[1] = ((pow(x[0], 2) + pow(y[0], 2)) * (x[2] - x[1]) 
           +(pow(x[1], 2) + pow(y[1], 2)) * (x[0] - x[2])
           +(pow(x[2], 2) + pow(y[2], 2)) * (x[1] - x[0])) / D;

    double A = 0.25 * sqrt((d[0] + d[1] + d[2]) * (-d[0] + d[1] + d[2]) * (d[0] - d[1] + d[2]) * (d[0] + d[1] - d[2]));
    *R = (d[0] * d[1] * d[2]) / (4.0 * A);

    // Derivatives

    // the gradients dRdv[:,i] are radial vectors, i.e. pointing alongside 
    //     dRdv[:,i] = (v[:,i] - c[:]) / |v[:,i] - c[:]| * lambda[i]
    // where lambda[i] is a scalar that we need to determine
    // we know that if we move all vertices radially outwarts by distance 1, then the circumradius increases by distance 1
    //     sum_i lambda[i] = 1
    // we do not have that lambda[i] > 0, i.e. moving a vertex towards the center can increase the circumradius (though this should only happen for obtuse triangles (?))

    // we now do some 2D geometry
    // we will move point A=0 radially outwards by distance eps and keep B=1, C=2, a=d[0] fixed
    // The circumcenter U will move along the perpendicular bisector of BC, by distance delta
    // Let M = midpoint(BC)
    // We have equations:
    //  |U - A| = |U - B| = |U - C| = R
    //  let h = |U - M|, then R^2 = |U - B|^2 = h^2 + (a/2)^2
    //  let theta = angle(MU, AU)
    //  the radius changes to R' = |U' - B| = sqrt((h+delta)^2 + (a/2)^2) = sqrt(R^2 + 2*h*delta) where delta^2=0 infinitesimally
    //  but also R' = |U' - A'| = |U-A + delta*(M-U)/h - eps*(A-U)/R|
    //  infinitesimally, we only care about the movement parallel to UA, which is of length (eps + delta*cos(theta))
    //  so we have dR/deps = 1 + delta/eps * cos(theta)
    //         and dR/deps = h/R * delta/eps
    //  so we have delta/eps = 1 / (h/R - cos(theta))
    //  and thus lambda = dR/deps = h/R / (h/R - cos(theta)) = 1 / (1 - cos(theta)*R/h)
    //  




    // // we now look at dcdv from the cartesian formula
    // // dc[0] / dv[0,0] = ((2*x[0]) * (y[1]-y[2])) / D - ((...)) / D^2 * dD/dv[0,0]
    // // where dD/dv[0,0] = 2*(y[1]-y[2])
    // // dc[1] / dv[0,0] = ((2*x[0]) * (x[2]-x[1]) + (x[1]²+y[1]²) * (+1) + (x[2]²+y[2]²) * (-1)) / D - ((...)) / D^2 * dD/dv[0,0]
    // // so we have
    // // dc[0] / dv[0,0] = 1/D * (2*x[0]*(y[1]-y[2]) - c[0]*2*(y[1]-y[2])) 
    // // = 2 * (y[1]-y[2]) * (x[0] - c[0]) / D 
    // // dc[1] / dv[0,0] = 1/D * (2*x[0]*(x[2]-x[1]) + (x[1]²+y[1]²) - (x[2]²+y[2]²) - c[1]*2*(y[1]-y[2]))
    // // = 1/D * (x[2]-x[1]) * (2*x[0] - x[2] - x[1]) + 1/D * (y[1]-y[2]) * (y[1]+y[2] - 2*c[1])

    // // for dc[*] / dv[1,0], we can use a formula that switches x and y everywhere; the sign of D, and thus the whole expression, also flips
    
    // for (int i = 0; i < 3; i++) {
    //     int j = (i + 1) % 3;
    //     int k = (i + 2) % 3;
    //     dcdv[0 * 6 + 0 * 3 + i] =   2.0 * (y[j] - y[k]) * (x[i] - c[0]) / D;
    //     dcdv[1 * 6 + 1 * 3 + i] = - 2.0 * (x[j] - x[k]) * (y[i] - c[1]) / D;
    //     dcdv[1 * 6 + 0 * 3 + i] =   (x[k] - x[j]) * (2.0 * x[i] - x[j] - x[k]) / D + (y[j] - y[k]) * (y[j] + y[k] - 2.0 * c[1]) / D;
    //     dcdv[0 * 6 + 1 * 3 + i] = - (y[k] - y[j]) * (2.0 * y[i] - y[j] - y[k]) / D + (x[j] - x[k]) * (x[j] + x[k] - 2.0 * c[0]) / D;
    // }




    // dcdv is not calculated yet
    return;
}

void welzl(
    // point set P
    const int nP,
    const double Px[], // [nP]
    const double Py[], // [nP]
    // boundary set R, in-out
    int *nR,
    double Rx[], // [3]
    double Ry[], // [3]
    // outputs
    double *cx, // center x
    double *cy, // center y
    double *R // radius
) {
    if (nR > 3) {
        error("welzl: nR > 3");
        return;
    } else if (nR == 3) { // trivial case, triangle
        // the smallest disk that contains 3 points is given either by
        // 1. the midpoint of the longest edge, with radius half the length of that edge; this is the case if the triangle is obtuse
        // 2. the circumcircle of the triangle formed by the 3 points

        // Edge lengths
        double d[3];
        d[0] = sqrt(pow(Rx[1] - Rx[2], 2) + pow(Ry[1] - Ry[2], 2));
        d[1] = sqrt(pow(Rx[0] - Rx[2], 2) + pow(Ry[0] - Ry[2], 2));
        d[2] = sqrt(pow(Rx[0] - Rx[1], 2) + pow(Ry[0] - Ry[1], 2));
        // Test for obtuse triangle
        for (int i = 0; i < 3; i++) {
            int j = (i + 1) % 3;
            int k = (i + 2) % 3;
            if (d[i] * d[i] > d[j] * d[j] + d[k] * d[k]) {
                // obtuse triangle, use midpoint of longest edge d[i], which is opposite vertex i
                *cx = 0.5 * (Rx[j] + Rx[k]);
                *cy = 0.5 * (Ry[j] + Ry[k]);
                *R = 0.5 * d[i];
                return;
            }
        }
        // acute or right triangle. points may be approximately degenerate so watch numerical stability
        // formula: https://en.wikipedia.org/wiki/Circumcircle

        // circumradius has a simple formula
        // R = (a*b*c) / (2*A)
        // where A is the area of the triangle
        // A = ((a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c))^0.5 / 4
        double A = 0.25 * sqrt((d[0] + d[1] + d[2]) * (-d[0] + d[1] + d[2]) * (d[0] - d[1] + d[2]) * (d[0] + d[1] - d[2]));
        *R = (d[0] * d[1] * d[2]) / (2.0 * A);

        // Barycentric coordinates
        //   U = a²(b² + c² - a²) : b²(c² + a² - b²) : c²(a² + b² - c²)
        // where a = |BC|, b = |AC|, c = |AB|
        // Using our index notation mod 3, the contribution of Rx[i] is
        //   d[i]²(d[i+1]² + d[i+2]² - d[i]²) = d[i]²( d2sum - 2*d[i]² )
        // The sum of the contributions is
        //   denominator = sum_i d[i]²( d2sum - 2*d[i]² ) = d2sum² - 2*sum_i d[i]^4
        double d2sum = pow(d[0], 2) + pow(d[1], 2) + pow(d[2], 2);
        double d4sum = pow(d[0], 4) + pow(d[1], 4) + pow(d[2], 4);
        double denominator = d2sum * d2sum - 2.0 * d4sum;
        // We now get the circumcenter
        double Ux = 0.0;
        double Uy = 0.0;
        for (int i = 0; i < 3; i++) {
            double coeff = pow(d[i], 2) * (d2sum - 2.0 * pow(d[i], 2)) / denominator;
            Ux += coeff * Rx[i];
            Uy += coeff * Ry[i];
        }
        *cx = Ux;
        *cy = Uy;
        return;
    } else if (nP == 0) { // 0, 1, or 2 boundary points
        if (nR == 0) { // no points => zero radius disk
            *cx = 0.0;
            *cy = 0.0;
            *R = 0.0;
            return;
        } else if (nR == 1) { // one point => zero radius disk at that point
            *cx = Rx[0];
            *cy = Ry[0];
            *R = 0.0;
            return;
        } else if (nR == 2) { // two points => disk at midpoint with radius half distance
            double mx = 0.5 * (Rx[0] + Rx[1]);
            double my = 0.5 * (Ry[0] + Ry[1]);
            double dx = Rx[0] - Rx[1];
            double dy = Ry[0] - Ry[1];
            double r = 0.5 * sqrt(dx * dx + dy * dy);
            *cx = mx;
            *cy = my;
            *R = r;
            return;
        }
    } else {
        // pick a random point from P. We just use the last one for simplicity
        double px = Px[nP - 1];
        double py = Py[nP - 1];
        // recurse without p
        welzl(nP - 1, Px, Py, nR, Rx, Ry, cx, cy, R);
        // if p is inside the circle, return
        double dx = px - *cx;
        double dy = py - *cy;
        if (dx * dx + dy * dy <= (*R) * (*R)) {
            return;
        }
        // otherwise, p is on the boundary
        int nR2 = nR + 1;
        double Rx2[3];
        double Ry2[3];
        memcpy(Rx2, Rx, nR * sizeof(double));
        memcpy(Ry2, Ry, nR * sizeof(double));
        Rx2[nR] = px;
        Ry2[nR] = py;
        welzl(nP - 1, Px, Py, nR2, Rx2, Ry2, cx, cy, R);
        return;
    }

    // Note: runtime is at worst 3 thrown away calls starting from nR=0 to nR=3
    // Each branch at worst recursed nP times
    // So we have 3*nP function calls, each of O(1) overhead => O(nP) total runtime
    // TODO: if this is a hotspot, we can remove the recursion for better CPU inline optimization
}

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