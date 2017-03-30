//// Body to field velocity influence coefficients
// Compute doublet and source velocity influence coefficient from one panel to a point
// Source induced velocities are given by Katz and Plotkin (program 12)
// Doublet induced velocities are computed from Katz and Plotkin (program 12) and smoothed with Vatistas vortex model
//
// Inputs:
// - wakeFlag: defines if wake or body panel
// - x1, y1, x2, y2, x3, y3, x4, y4: coordinates of infulencing panel vertices in panel axis
// - x, y, z: coordinates of target point in panel axis
// - l0, l1, l2: components of influencing panel longitudinal unit vector
// - p0, p1, p2: components of influencing panel perpendicular unit vector
// - n0, n1, n2: components of influencing panel normal unit vector
//
// Output:
// - coeff: array of AIC ([0] = mu-u, [1] = mu-v, [2] = mu-w, [3] = tau-u, [4] = tau-v, [5] = tau-w)

#include <iostream>
#include <array>
#include <Eigen/Dense>
#include "infcBF.h"

#define PI 3.14159
#define ALPHA 0.2 // ALPHA is the smoothing parameter for doublets. Its value should be ~[.2, .1] (check refs)
#define BETA 0.05 // BETA is the smoothing parameter for source (U,V). Its value should be ~[.05, .1] (check refs)

using namespace std;
using namespace Eigen;

array<double,6> infcBF(bool wakeFlag, double x1 ,double x2, double x3, double x4, double y1, double y2, double y3, double y4,
                       double x, double y, double z, double l0, double l1, double l2,
                       double p0, double p1, double p2, double n0, double n1, double n2) {

    // Temporary variables
    RowVector4d xV, yV; // Panel vertices (clockwise order)
    RowVector4d d, r, e, h, F, G; // Katz and Plotkin coefficients
    RowVector4d d1, d2, r1_, r2_, h_c;
    MatrixXd r1V, r2V; // Vertex to point vectors for a panel edge
    RowVector4d hV, r_c;
    r1V.resize(3,4);
    r2V.resize(3,4);

    double muU = 0; // ICs
    double muV = 0;
    double muW = 0;
    double tauU = 0;
    double tauV = 0;
    double tauW = 0;
    RowVector4d dummy;
    dummy(0) = 0;
    dummy(1) = 0;
    dummy(2) = 0;
    dummy(3) = 0;
    array<double,6> coeff;

    // Local copy of coordinates
    xV(0) = x1; xV(1) = x2; xV(2) = x3; xV(3) = x4;
    yV(0) = y1; yV(1) = y2; yV(2) = y3; yV(3) = y4;

    //// Coefficients
    // Genral
    for (int i = 1; i <= 4; ++i) {
        int j = i%4 + 1;

        d(i-1) = sqrt((xV(j-1) - xV(i-1))*(xV(j-1) - xV(i-1)) + (yV(j-1) - yV(i-1))*(yV(j-1) - yV(i-1)));
        r(i-1) = sqrt((x - xV(i-1))*(x - xV(i-1)) + (y - yV(i-1))*(y - yV(i-1)) + z * z);
        e(i-1) = (x - xV(i-1))*(x - xV(i-1)) + z*z;
        h(i-1) = (x - xV(i-1)) * (y - yV(i-1));
    }
    for (int i = 1; i <= 4; ++i) {
        int j = i%4 + 1;

        F(i-1) = (yV(j-1) - yV(i-1))*e(i-1) - (xV(j-1) - xV(i-1))*h(i-1);
        G(i-1) = (yV(j-1) - yV(i-1))*e(j-1) - (xV(j-1) - xV(i-1))*h(j-1);
    }
    // Smoothing
    for (int i = 1; i <= 4; ++i) {
        int j = i % 4 + 1;

        r1V(0,i-1) = x - xV(i-1); r1V(1,i-1) = y - yV(i-1); r1V(2,i-1) = z;
        r2V(0,i-1) = x - xV(j-1); r2V(1,i-1) = y - yV(j-1); r2V(2,i-1) = z;
        hV(i-1) = sqrt((r1V(1,i-1)*r2V(2,i-1)-r1V(2,i-1)*r2V(1,i-1))*(r1V(1,i-1)*r2V(2,i-1)-r1V(2,i-1)*r2V(1,i-1))
                     + (r1V(2,i-1)*r2V(0,i-1)-r1V(0,i-1)*r2V(2,i-1))*(r1V(2,i-1)*r2V(0,i-1)-r1V(0,i-1)*r2V(2,i-1))
                     + (r1V(0,i-1)*r2V(1,i-1)-r1V(1,i-1)*r2V(0,i-1))*(r1V(0,i-1)*r2V(1,i-1)-r1V(1,i-1)*r2V(0,i-1)));
        hV(i-1) = hV(i-1) / d(i-1);
        r_c(i-1) = ALPHA * d(i-1);
        }
    for (int i = 1; i <= 4; ++i) {
        int j = i % 4 + 1;

        h_c(i-1) = BETA * d(i-1);
        d1(i-1) = r(i-1) * r(i-1) - hV(i-1) * hV(i-1);
        d2(i-1) = r(j-1) * r(j-1) - hV(i-1) * hV(i-1);
        r1_(i-1) = sqrt(h_c(i-1) * h_c(i-1) + d1(i-1));
        r2_(i-1) = sqrt(h_c(i-1) * h_c(i-1) + d2(i-1));
    }

    //// AIC
    // Source (from Katz & Plotkin, program 12)
    if (!wakeFlag) {
        for (int i = 1; i <= 4; ++i) {
            int j = i % 4 + 1;

            if (hV(i-1) < h_c(i-1)) {
                if (i == 2 || i == 4)
                    tauU = tauU + (yV(j - 1) - yV(i - 1)) / d(i - 1) *
                                  log((r1_(i - 1) + r2_(i - 1) - d(i - 1)) / (r1_(i - 1) + r2_(i - 1) + d(i - 1)));
                else
                    tauU = tauU + (yV(j - 1) - yV(i - 1)) / d(i - 1) *
                                  log((r(i - 1) + r(j - 1) - d(i - 1)) / (r(i - 1) + r(j - 1) + d(i - 1)));
                if (i == 1 || i == 3)
                    tauV = tauV + (xV(j - 1) - xV(i - 1)) / d(i - 1) *
                                  log((r1_(i - 1) + r2_(i - 1) - d(i - 1)) / (r1_(i - 1) + r2_(i - 1) + d(i - 1)));
                else
                    tauV = tauV + (xV(j - 1) - xV(i - 1)) / d(i - 1) *
                                  log((r(i - 1) + r(j - 1) - d(i - 1)) / (r(i - 1) + r(j - 1) + d(i - 1)));

                tauW = tauW + atan2(z * (xV(j-1) - xV(i-1)) * (F(i-1) * r(j-1) - G(i-1) * r(i-1)),
                                    z * z * (xV(j-1) - xV(i-1)) * (xV(j-1) - xV(i-1)) * r(i-1) * r(j-1) + F(i-1) * G(i-1));
            }
            else {
                tauU = tauU + (yV(j - 1) - yV(i - 1)) / d(i - 1) *
                              log((r(i - 1) + r(j - 1) - d(i - 1)) / (r(i - 1) + r(j - 1) + d(i - 1)));
                tauV = tauV + (xV(j - 1) - xV(i - 1)) / d(i - 1) *
                              log((r(i - 1) + r(j - 1) - d(i - 1)) / (r(i - 1) + r(j - 1) + d(i - 1)));
                tauW = tauW + atan2(z * (xV(j-1) - xV(i-1)) * (F(i-1) * r(j-1) - G(i-1) * r(i-1)),
                                    z * z * (xV(j-1) - xV(i-1)) * (xV(j-1) - xV(i-1)) * r(i-1) * r(j-1) + F(i-1) * G(i-1));
            }
        }
        tauU = 1/(4*PI) * tauU;
        tauV = -1/(4*PI) * tauV;
        tauW = 1/(4*PI) * tauW;
    }
    else {
        tauU = 0;
        tauV = 0;
        tauW = 0;
    }

    // Doublet (from Katz & Plotkin + Vatistas smoothing, n = 1)
    for (int i = 1; i <= 4; ++i) {
        int j = i % 4 + 1;

        dummy(0) = (y - yV(i-1)) * z - z * (y - yV(j-1));
        dummy(1) = z * (x - xV(j-1)) - (x - xV(i-1)) * z;
        dummy(2) = (x - xV(i-1)) * (y - yV(j-1)) - (y - yV(i-1)) * (x - xV(j-1));
        dummy(3) = (r(i-1)+r(j-1))
                   / (r(i-1)*r(j-1)*(r(i-1)*r(j-1) + ((x-xV(i-1))*(x-xV(j-1)) + (y-yV(i-1))*(y-yV(j-1)) + z*z)));

        muU = muU + hV(i-1)*hV(i-1) / (r_c(i-1)*r_c(i-1) + hV(i-1)*hV(i-1)) * dummy(0) * dummy(3);
        muV = muV + hV(i-1)*hV(i-1) / (r_c(i-1)*r_c(i-1) + hV(i-1)*hV(i-1)) * dummy(1) * dummy(3);
        muW = muW + hV(i-1)*hV(i-1) / (r_c(i-1)*r_c(i-1) + hV(i-1)*hV(i-1)) * dummy(2) * dummy(3);
    }
    muU = -1/(4*PI) * muU;
    muV = -1/(4*PI) * muV;
    muW = -1/(4*PI) * muW;

    coeff[0] =  muU * l0 + muV * p0 + muW * n0;
    coeff[1] =  muU * l1 + muV * p1 + muW * n1;
    coeff[2] =  muU * l2 + muV * p2 + muW * n2;
    coeff[3] =  tauU * l0 + tauV * p0 + tauW * n0;
    coeff[4] =  tauU * l1 + tauV * p1 + tauW * n1;
    coeff[5] =  tauU * l2 + tauV * p2 + tauW * n2;
    return coeff;
}