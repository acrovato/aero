//// Split panel
// Split a panel and compute AIC on each subpanel
//

#include <iostream>
#include <Eigen/Dense>
#include <array>
#include "split_panel.h"
#include "infcBF.h"

using namespace std;
using namespace Eigen;

array<RowVectorXd,6> split_panel(int NS, int NSs, double x1, double x2, double x3, double x4,
                                 double y1, double y2, double y3, double y4,
                                 double x, double y, double z,
                                 double l0, double l1, double l2,
                                 double p0, double p1, double p2,
                                 double n0, double n1, double n2) {

    // Initialization
    int idx = 0; // counter
    double a0, a1 ,b0, b1; // interpolation coefficients
    double xP1, xP2, xP3, xP4, yP1, yP2, yP3, yP4; // sub-panel vertices
    array <RowVectorXd, 6> coeff; // returned and temporary coefficients
    array <double, 6> coeffT;
    coeff[0].resize(NS);
    coeff[1].resize(NS);
    coeff[2].resize(NS);
    coeff[3].resize(NS);
    coeff[4].resize(NS);
    coeff[5].resize(NS);

    // Split panel
    // TODO: consider storing xP, yP to find closest subpanel to field point and shifting field point above that subpanel
    for (int j = 0; j < NSs; j++) {
        for (int i = 0; i < NSs; i++) {
            // Compute weight factors
            a0 = (double) i/NSs;
            a1 = (double) (i+1)/NSs;
            b0 = (double) j/NSs;
            b1 = (double) (j+1)/NSs;
            // Compute new corner points
            xP1 = (1-b0)*((1-a0)*x1 + a0*x2) + b0*(a0*x3 +(1-a0)*x4);
            xP2 = (1-b0)*((1-a1)*x1 + a1*x2) + b0*(a1*x3 +(1-a1)*x4);
            xP3 = (1-b1)*((1-a1)*x1 + a1*x2) + b1*(a1*x3 +(1-a1)*x4);
            xP4 = (1-b1)*((1-a0)*x1 + a0*x2) + b1*(a0*x3 +(1-a0)*x4);
            yP1 = (1-b0)*((1-a0)*y1 + a0*y2) + b0*(a0*y3 +(1-a0)*y4);
            yP2 = (1-b0)*((1-a1)*y1 + a1*y2) + b0*(a1*y3 +(1-a1)*y4);
            yP3 = (1-b1)*((1-a1)*y1 + a1*y2) + b1*(a1*y3 +(1-a1)*y4);
            yP4 = (1-b1)*((1-a0)*y1 + a0*y2) + b1*(a0*y3 +(1-a0)*y4);
            // call infcBF
            coeffT = infcBF(0, xP1, xP2, xP3, xP4, yP1, yP2, yP3, yP4,
                            x, y, z, l0, l1, l2, p0, p1, p2, n0, n1, n2);
            // Store to return
            coeff[0](idx) = coeffT[0];
            coeff[1](idx) = coeffT[1];
            coeff[2](idx) = coeffT[2];
            coeff[3](idx) = coeffT[3];
            coeff[4](idx) = coeffT[4];
            coeff[5](idx) = coeffT[5];
            idx ++;
        }
    }
    return coeff;
}
