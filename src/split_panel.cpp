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
#include "infcBF.h"

#define NDIM 3
#define NMG 7
#define NSING 6

using namespace std;
using namespace Eigen;

array<RowVectorXd,42> split_panel(int idP, int idF,
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
    MatrixX3d tgt; // target points
    tgt.resize(NMG, NDIM);

    array <RowVectorXd, 42> coeff; // returned and temporary coefficients
    array <double, NSING> coeffT;
    for (int i = 0; i < 42; i++)
        coeff[i].resize(sp.NS);

    //// Change coordinates
    // Cell center
    tgt(0,0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(0,1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(0,2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0)) + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));
    // X-bwd
    tgt(1,0) = bPan.l(idP,0)*(fPan.CG(idF,0)-fPan.deltaMG) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(1,1) = bPan.p(idP,0)*(fPan.CG(idF,0)-fPan.deltaMG) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(1,2) = bPan.n(idP,0)*(fPan.CG(idF,0)-fPan.deltaMG - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));
    // X-fwd
    tgt(2,0) = bPan.l(idP,0)*(fPan.CG(idF,0)+fPan.deltaMG) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(2,1) = bPan.p(idP,0)*(fPan.CG(idF,0)+fPan.deltaMG) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(2,2) = bPan.n(idP,0)*(fPan.CG(idF,0)+fPan.deltaMG - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));
    // Y-bwd
    tgt(3,0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*(fPan.CG(idF,1)-fPan.deltaMG) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(3,1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*(fPan.CG(idF,1)-fPan.deltaMG) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(3,2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1)-fPan.deltaMG - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));
    // Y-fwd
    tgt(4,0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*(fPan.CG(idF,1)+fPan.deltaMG) + bPan.l(idP,2)*fPan.CG(idF,2);
    tgt(4,1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*(fPan.CG(idF,1)+fPan.deltaMG) + bPan.p(idP,2)*fPan.CG(idF,2);
    tgt(4,2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1)+fPan.deltaMG - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2) - bPan.CG(idP,2));
    // Z-bwd
    tgt(5,0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*(fPan.CG(idF,2)-fPan.deltaMG);
    tgt(5,1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*(fPan.CG(idF,2)-fPan.deltaMG);
    tgt(5,2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2)-fPan.deltaMG - bPan.CG(idP,2));
    // Z-fwd
    tgt(6,0) = bPan.l(idP,0)*fPan.CG(idF,0) + bPan.l(idP,1)*fPan.CG(idF,1) + bPan.l(idP,2)*(fPan.CG(idF,2)+fPan.deltaMG);
    tgt(6,1) = bPan.p(idP,0)*fPan.CG(idF,0) + bPan.p(idP,1)*fPan.CG(idF,1) + bPan.p(idP,2)*(fPan.CG(idF,2)+fPan.deltaMG);
    tgt(6,2) = bPan.n(idP,0)*(fPan.CG(idF,0) - bPan.CG(idP,0))
               + bPan.n(idP,1)*(fPan.CG(idF,1) - bPan.CG(idP,1))
               + bPan.n(idP,2)*(fPan.CG(idF,2)+fPan.deltaMG - bPan.CG(idP,2));

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
    for (int i = 0; i < NMG; i++) {
        for (int j = 0; j < sp.NS; j++) {
            // call infcBB
            coeffT = infcBF(0, xV0(j), xV1(j), xV2(j), xV3(j), yV0(j), yV1(j), yV2(j), yV3(j),
                            tgt(i,0), tgt(i,1), tgt(i,2), bPan.l(idP,0), bPan.l(idP,1), bPan.l(idP,2),
                            bPan.p(idP,0), bPan.p(idP,1), bPan.p(idP,2), bPan.n(idP,0), bPan.n(idP,1), bPan.n(idP,2));
            // Store to return
            coeff[i*NSING+0](j) = coeffT[0];
            coeff[i*NSING+1](j) = coeffT[1];
            coeff[i*NSING+2](j) = coeffT[2];
            coeff[i*NSING+3](j) = coeffT[3];
            coeff[i*NSING+4](j) = coeffT[4];
            coeff[i*NSING+5](j) = coeffT[5];
        }
    }
    return coeff;
}
