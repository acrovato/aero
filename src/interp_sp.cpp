//// Interpolation from panel to sub-panel
// Interpolate linearly surface singularities from panel vertices to sub-panels
//
// Inputs:
// - idP: current panel index
// - bPan: body panels (structure)
// - mu0, mu1, mu2, mu3: doublet at interpolation points (vertices)
// - tau0, tau1, tau2, tau3: sources at interpolation points (vertices)
//
// Output:
// - spis: interpolated singularities (row = sub-panel number; col 0 = doublet, col 1 = source)

#include <iostream>
#include <Eigen/Dense>
#include "interp_sp.h"
#include "interp.h"

using namespace Eigen;

MatrixXd interp_sp(int idP, Network &bPan, Subpanel &sp,
                   double mu0, double mu1, double mu2, double mu3,
                   double tau0, double tau1, double tau2, double tau3) {

    // Temporary variables
    double x0, x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3; // vertices coordinates
    double a, b; // interpolation parameters
    double X, Y, Z; // sub-panel center
    MatrixXd spis; // Interpolated singularities
    spis.resize(sp.NS, 2);

    // Copy vertices ton local variables
    x0 = bPan.v0(idP,0);
    x1 = bPan.v1(idP,0);
    x2 = bPan.v2(idP,0);
    x3 = bPan.v3(idP,0);
    y0 = bPan.v0(idP,1);
    y1 = bPan.v1(idP,1);
    y2 = bPan.v2(idP,1);
    y3 = bPan.v3(idP,1);
    z0 = bPan.v0(idP,2);
    z1 = bPan.v1(idP,2);
    z2 = bPan.v2(idP,2);
    z3 = bPan.v3(idP,2);

    // Computation of sub-panel centers in global axes
    int idx = 0, j = 0;
    for (int jj = 0; jj < sp.NSs; jj++) {
        int i = 0;
        for (int ii = 0; ii < sp.NSs; ii++) {
            // Compute weight factors
            a = (double) (i+1)/sp.NSs/2;
            b = (double) (j+1)/sp.NSs/2;
            // Compute center point
            X = (1-b)*((1-a)*x0 + a*x1) + b*(a*x2 + (1-a)*x3);
            Y = (1-b)*((1-a)*y0 + a*y1) + b*(a*y2 + (1-a)*y3);
            Z = (1-b)*((1-a)*z0 + a*z1) + b*(a*z2 + (1-a)*z3);

            // Interpolate singularities
            spis(idx,0) = interp(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3,
                                 mu0, mu1, mu2, mu3,
                                 X, Y, Z);
            spis(idx,1) = interp(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3,
                                 tau0, tau1, tau2, tau3,
                                 X, Y, Z);
            idx++;
            i += 2;
        }
        j += 2;
    }
    return spis;
}