//// Interpolation
// Interpolate linearly a quantity defined in 4 points to another point
//
// Inputs:
// - x0, x1, x2, x3: x-coordinates of interpolation points
// - y0, y1, y2, y3: y-coordinates of interpolation points
// - z0, z1, z2, z3: z-coordinates of interpolation points
// - s0, s1, s2, s3: quantity to interpolate at interpolation point
// - x, y, z: coordinates of interpolated point
//
// Output:
// - i: interpolated quantity

#include <math.h>
#include "interp.h"

double interp(double x0, double y0, double z0, double x1, double y1, double z1,
              double x2, double y2, double z2, double x3, double y3, double z3,
              double s0, double s1, double s2, double s3,
              double x, double y, double z) {

    // Initialization
    double u0, u1, t;
    double i0, i1, i;

    // Compute modified distance ratios for interpolation
    u0 = sqrt((x-x0)*(x-x0) + (z-z0)*(z-z0)) / sqrt((x1-x0)*(x1-x0) + (z1-z0)*(z1-z0));
    u1 = sqrt((x-x3)*(x-x3) + (z-z3)*(z-z3)) / sqrt((x2-x3)*(x2-x3) + (z2-z3)*(z2-z3));
    t = sqrt((y-y0)*(y-y0) + (z-z0)*(z-z0)) / sqrt((y3-y0)*(y3-y0) + (z3-z0)*(z3-z0));

    // Interpolate along x
    i0 = (1-u0)*s0 + u0*s1;
    i1 = u1*s2 + (1-u1)*s3;

    // Interpolate along y
    i = (1-t)*i0 + t*i1;
    return i;
}