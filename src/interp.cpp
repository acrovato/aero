//// Bilinear interpolation
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

#include <iostream>
#include <Eigen/Dense>
#include "interp.h"

#define TOL 1e-2 // geometric tolerance on z (measure how far target point is w.r.t. interpolation plane)

using namespace std;
using namespace Eigen;

double interp(double x0, double y0, double z0, double x1, double y1, double z1,
              double x2, double y2, double z2, double x3, double y3, double z3,
              double s0, double s1, double s2, double s3,
              double x, double y, double z) {

    //// Initialization
    Vector3d l, p, n; // unit vector
    Vector3d cg; // CG of planed formed by interpolation points
    Vector2d A, B, C, D, X; // (rotated) interpolation and interpolated points
    Vector2d E, F, G, H; // temporary vector
    double zd, avg; // distance between interpolated point and interpolation plane & plane diagonal
    double k2, k1, k0, Delta; // coefficients of k2*x^2 + k1*x + k0 and discriminant
    double u, v ;// inverse bilinear interpolation parameters
    double i0, i1, i ;// interpolated quantities

    //// TODO if it works, move change of coordinates out of file to improve CPU runtime
    //// Change of axes
    // Plane center
    cg(0) = 0.25*(x0 + x1 + x2 + x3);
    cg(1) = 0.25*(y0 + y1 + y2 + y3);
    cg(2) = 0.25*(z0 + z1 + z2 + z3);
    // Unit vectors
    l(0) = 0.5*(0.5*(x0+x3)-0.5*(x1+x2));
    l(1) = 0.5*(0.5*(y0+y3)-0.5*(y1+y2));
    l(2) = 0.5*(0.5*(z0+z3)-0.5*(z1+z2));
    n(0) = (y2-y0)*(z1-z3) - (z2-z0)*(y1-y3);
    n(1) = (z2-z0)*(x1-x3) - (x2-x0)*(z1-z3);
    n(2) = (x2-x0)*(y1-y3) - (y2-y0)*(x1-x3);
    p = n.cross(l);
    l = l / l.norm(); // normalize
    p = p / p.norm(); // normalize
    n = n / n.norm(); // normalize
    // Translation & Rotation
    A(0) = l(0)*(cg(0)-x0) + l(1)*(cg(1)-y0) + l(2)*(cg(2)-z0);
    A(1) = p(0)*(cg(0)-x0) + p(1)*(cg(1)-y0) + p(2)*(cg(2)-z0);
    B(0) = l(0)*(cg(0)-x1) + l(1)*(cg(1)-y1) + l(2)*(cg(2)-z1);
    B(1) = p(0)*(cg(0)-x1) + p(1)*(cg(1)-y1) + p(2)*(cg(2)-z1);
    C(0) = l(0)*(cg(0)-x2) + l(1)*(cg(1)-y2) + l(2)*(cg(2)-z2);
    C(1) = p(0)*(cg(0)-x2) + p(1)*(cg(1)-y2) + p(2)*(cg(2)-z2);
    D(0) = l(0)*(cg(0)-x3) + l(1)*(cg(1)-y3) + l(2)*(cg(2)-z3);
    D(1) = p(0)*(cg(0)-x3) + p(1)*(cg(1)-y3) + p(2)*(cg(2)-z3);
    X(0) = l(0)*(cg(0)-x) + l(1)*(cg(1)-y) + l(2)*(cg(2)-z);
    X(1) = p(0)*(cg(0)-x) + p(1)*(cg(1)-y) + p(2)*(cg(2)-z);
    // Error measure
    avg = sqrt((x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) + (z2-z0)*(z2-z0));
    zd = n(0)*(x-cg(0)) + n(1)*(y-cg(1)) + n(2)*(z-cg(2));
    if (abs(zd/avg) > TOL) {
        cout << endl << "Interpolated point is far from interpolation plane; bilinear interpolation might be inaccurate!" << endl;
        cout << "Distance versus interpolation plane diagonal: " << zd << " / " << avg << " > " << TOL << endl;
        cout << "This might be caused by a coarse panelling of the body surface." << endl;
    }

    //// Interpolation
    // Intermediate vectors
    E = B-A;
    F = D-A;
    G = (A-B) + (C-D);
    H = X-A;
    // Coefficients
    k2 = G(0)*F(1) - G(1)*F(0); //
    k1 = E(0)*F(1) - E(1)*F(0) + H(0)*G(1) - H(1)*G(0);
    k0 = H(0)*E(1) - H(1)*E(0);
    // Solve quadratic equation
    if (k2 == 0) {
        v = -k0/k1;
    }
    else {
        Delta = k1*k1 - 4*k2*k0;
        if (Delta < 0) {
            cout << endl << "Discriminant of inverse bilinear interpolation is negative. Cannot interpolate!" << endl;
            exit(EXIT_FAILURE);
        }
        else if (Delta == 0) {
            v = (-k1) / (2.0*k2);
        }
        else {
            Delta = sqrt(Delta);
            v = (-k1 - Delta) / (2.0*k2);
        }
    }
    u = (H(0) - F(0)*v) / (E(0) + G(0)*v);
    // Interpolation along l
    i0 = (1-u)*s0 + u*s1;
    i1 = (1-u)*s3 + u*s2;
    // Interpolation along p
    i = (1-v)*i0 + v*i1;
    return i;
}