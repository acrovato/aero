//// Body to body/field potential influence coefficient computation
// Compute doublet and source potential influence coefficient between 2 surface panels
// Katz and Plotkin (2001), Low-Speed Aerodynamics. Program 12
//
// Inputs:
// - wakeFlag: defines if wake or body panel
// - idTgt: index of target panel
// - idSrc: index of influencing panel
// - x, y, z: coordinates of target panel center in panel axis
// - x1, y1, x2, y2, x3, y3, x4, y4: coordinates of influencing panel vertices in panel axis
//
// Output:
// - coeff: array of AIC ([0] = mu, [1] = tau)

#include <iostream>
#include <Eigen/Dense>
#include <array>
#include "infcB.h"

#define PI 3.14159

using namespace std;
using namespace Eigen;

array<double,2> infcB(bool wakeFlag,int idTgt, int idSrc, double x, double y, double z,
                      double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {

    // Temporary
    RowVector4d xV, yV; // Local vertices
    RowVector4d d, r, e, h, F, G; // coefficients
    double mu = 0, tau = 0; // AIC
    array<double,2> coeff;

    xV(0) = x1; xV(1) = x2; xV(2) = x3; xV(3) = x4;
    yV(0) = y1; yV(1) = y2; yV(2) = y3; yV(3) = y4;

    //// Wake panels
    if (wakeFlag) {

        // Coefficients
        for (int i = 1; i <= 4; ++i) {
            int j = i%4 + 1;

            r(i-1) = sqrt((x - xV(i-1))*(x - xV(i-1)) + (y - yV(i-1))*(y - yV(i-1)) + z * z);
            e(i-1) = (x - xV(i-1))*(x - xV(i-1)) + z*z;
            h(i-1) = (x - xV(i-1)) * (y - yV(i-1));
        }
        for (int i = 1; i <= 4; ++i) {
            int j = i%4 + 1;

            F(i-1) = (yV(j-1) - yV(i-1))*e(i-1) - (xV(j-1) - xV(i-1))*h(i-1);
            G(i-1) = (yV(j-1) - yV(i-1))*e(j-1) - (xV(j-1) - xV(i-1))*h(j-1);
        }

        // Doublet
        for (int i = 1; i <= 4; ++i) {
            int j = i % 4 + 1;

            mu = mu + atan2(z * (xV(j-1) - xV(i-1)) * (F(i-1) * r(j-1) - G(i-1) * r(i-1)),
                            z * z * (xV(j-1) - xV(i-1)) * (xV(j-1) - xV(i-1)) * r(i-1) * r(j-1) + F(i-1) * G(i-1));
        }
        mu = -1 / (4 * PI) * mu;
        tau = 0;
    }

    //// Body panels
    else {
        // Coefficients
        for (int i = 1; i <= 4; ++i) {
            int j = i%4 + 1;

            d(i-1) = sqrt((xV(j-1) - xV(i-1))*(xV(j-1) - xV(i-1)) + (yV(j-1) - yV(i-1))*(yV(j-1) - yV(i-1)));
            r(i-1) = sqrt((x - xV(i-1))*(x - xV(i-1)) + (y - yV(i-1))*(y - yV(i-1)) + z * z);
        }

        if (idTgt == idSrc) {
            // Source
            for (int i = 1; i <= 4; ++i) {
                int j = i % 4 + 1;

                tau = tau + ((x - xV(i-1)) * (yV(j-1) - yV(i-1)) - (y - yV(i-1)) * (xV(j-1) - xV(i-1))) / d(i-1) *
                            log((r(i-1) + r(j-1) + d(i-1)) / (r(i-1) + r(j-1) - d(i-1)));
            }
            mu = 0.5;
            tau = -1 / (4 * PI) * tau;
        }
        else {
            // Coefficients
            for (int i = 1; i <= 4; ++i) {
                e(i-1) = (x - xV(i-1))*(x - xV(i-1)) + z*z;
                h(i-1) = (x - xV(i-1)) * (y - yV(i-1));
            }
            for (int i = 1; i <= 4; ++i) {
                int j = i%4 + 1;

                F(i-1) = (yV(j-1) - yV(i-1))*e(i-1) - (xV(j-1) - xV(i-1))*h(i-1);
                G(i-1) = (yV(j-1) - yV(i-1))*e(j-1) - (xV(j-1) - xV(i-1))*h(j-1);
            }

            // Doublet
            for (int i = 1; i <= 4; ++i) {
                int j = i % 4 + 1;

                mu = mu + atan2(z * (xV(j-1) - xV(i-1)) * (F(i-1) * r(j-1) - G(i-1) * r(i-1)),
                                      z * z * (xV(j-1) - xV(i-1)) * (xV(j-1) - xV(i-1)) * r(i-1) * r(j-1) + F(i-1) * G(i-1));
            }
            mu = -1 / (4 * PI) * mu;
            // Source
            for (int i = 1; i <= 4; ++i) {
                int j = i % 4 + 1;

                tau = tau + ((x - xV(i-1)) * (yV(j-1) - yV(i-1)) - (y - yV(i-1)) * (xV(j-1) - xV(i-1))) / d(i-1) *
                                log((r(i-1) + r(j-1) + d(i-1)) / (r(i-1) + r(j-1) - d(i-1)));
            }
            tau = -1 / (4 * PI) * tau - z * mu;
        }
    }
    coeff[0] = mu;
    coeff[1] = tau;

    return  coeff;
}