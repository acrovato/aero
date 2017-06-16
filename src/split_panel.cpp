//// Split panel
// Split a panel and compute AIC on each subpanel
//
// Inputs:
// - idP: influencing panel index
// - idF: target field cell index
// - x0: x - coordinate of first corner point of panel in local axes
// - x1: x - coordinate of second corner point of panel in local axes
// - x2: x - coordinate of third corner point of panel in local axes
// - x3: x - coordinate of fourth corner point of panel in local axes
// - y0: y - coordinate of first corner point of panel in local axes
// - y1: y - coordinate of second corner point of panel in local axes
// - y2: y - coordinate of third corner point of panel in local axes
// - y3: y - coordinate of fourth corner point of panel in local axes
// - bPan: body panels (structure)
// - fPan: field panels (structure)
// - sp: sub-panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include <array>
#include "split_panel.h"
#include "infcB.h"

#define NDIM 3
#define NSING 2

using namespace std;
using namespace Eigen;

array<RowVectorXd,NSING> split_panel(int idP, int idF,
                                  double x0, double x1, double x2, double x3, double y0, double y1, double y2, double y3,
                                  Network &bPan, Field &fPan, Subpanel &sp) {

    //// Initialization
    int idx; // counter
    double a0, a1 ,b0, b1; // interpolation coefficients
    VectorXd xV0, xV1, xV2, xV3, yV0, yV1, yV2, yV3; // sub-panel vertices
    xV0.resize(sp.NS);
    xV1.resize(sp.NS);
    xV2.resize(sp.NS);
    xV3.resize(sp.NS);
    yV0.resize(sp.NS);
    yV1.resize(sp.NS);
    yV2.resize(sp.NS);
    yV3.resize(sp.NS);
    Vector3d tgt; // target points
    tgt.resize(NDIM);

    array <RowVectorXd, NSING> coeff; // returned and temporary coefficients
    array <double, NSING> coeffT;
    for (int i = 0; i < NSING; i++)
        coeff[i].resize(sp.NS);

    //// Change coordinates
    // Cell center
    tgt(0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0)) + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));

    //// Split panel
    // Compute sub-panels vertices
    idx = 0;
    for (int j = 0; j < sp.NSs; j++) {
        for (int i = 0; i < sp.NSs; i++) {
            // Compute weight factors
            a0 = (double) i/sp.NSs;
            a1 = (double) (i+1)/sp.NSs;
            b0 = (double) j/sp.NSs;
            b1 = (double) (j+1)/sp.NSs;
            // Compute new vertices
            xV0(idx) = (1-b0)*((1-a0)*x0 + a0*x1) + b0*(a0*x2 +(1-a0)*x3);
            xV1(idx) = (1-b0)*((1-a1)*x0 + a1*x1) + b0*(a1*x2 +(1-a1)*x3);
            xV2(idx) = (1-b1)*((1-a1)*x0 + a1*x1) + b1*(a1*x2 +(1-a1)*x3);
            xV3(idx) = (1-b1)*((1-a0)*x0 + a0*x1) + b1*(a0*x2 +(1-a0)*x3);
            yV0(idx) = (1-b0)*((1-a0)*y0 + a0*y1) + b0*(a0*y2 +(1-a0)*y3);
            yV1(idx) = (1-b0)*((1-a1)*y0 + a1*y1) + b0*(a1*y2 +(1-a1)*y3);
            yV2(idx) = (1-b1)*((1-a1)*y0 + a1*y1) + b1*(a1*y2 +(1-a1)*y3);
            yV3(idx) = (1-b1)*((1-a0)*y0 + a0*y1) + b1*(a0*y2 +(1-a0)*y3);
            idx ++;
        }
    }

    //// Compute AIC
    for (int j = 0; j < sp.NS; j++) {
        // call infcBB
        coeffT = infcB(0, 0, 1, tgt(0), tgt(1), tgt(2), xV0(j), yV0(j), xV1(j), yV1(j), xV2(j), yV2(j), xV3(j), yV3(j));
        // Store to return
        coeff[0](j) = coeffT[0];
        coeff[1](j) = coeffT[1];
    }
    return coeff;
}
