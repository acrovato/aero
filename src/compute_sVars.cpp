// Flow Computation
// Compute flow physical properties (body velocity and pressure and field mach)
// Currently computes body velocity by FD, incompressible pressure coefficient and forces

#include <iostream>
#include <Eigen/Dense>

#include "compute_sVars.h"

#define NDIM 3
#define GAMMA 1.4

using namespace std;
using namespace Eigen;

void compute_sVars(bool symY, double sRef, double alpha, double machInf, Vector3d &vInf,
               MatrixX3d &vSigma, Network &bPan, double &cL, double &cD) {

    // Temporary variables
    int panIdx = 0;
    double deltaL1, deltaL2, deltaT1, deltaT2;
    double u = 0, v = 0, w = 0;
    double cX = 0, cZ = 0;

    //// Begin
    cout << "Computing surface variables... ";

    // Resize matrices
    bPan.U = MatrixX3d::Zero(bPan.nP, NDIM);
    bPan.cP.resize(bPan.nP, 1);

    //// Velocity computation
    // Surface induced velocity
    for (int j = 0; j < bPan.nS_; ++j) {
        for (int i = 0; i < bPan.nC_; ++i) {

            // Longitudinal
            if (i == 0) { // Forward differencing
                deltaL2 = (bPan.CG.row(panIdx + 1) - bPan.CG.row(panIdx)).norm();
                u = (bPan.mu(panIdx + 1) - bPan.mu(panIdx)) / deltaL2;
            }
            else if (i == bPan.nC_ - 1) { // Backward differencing
                deltaL1 = (bPan.CG.row(panIdx) - bPan.CG.row(panIdx - 1)).norm();
                u = (bPan.mu(panIdx) - bPan.mu(panIdx - 1)) / deltaL1;
            }
            else { // Central differencing
                deltaL1 = (bPan.CG.row(panIdx) - bPan.CG.row(panIdx - 1)).norm();
                deltaL2 = (bPan.CG.row(panIdx + 1) - bPan.CG.row(panIdx)).norm();
                u = 0.5 * ((bPan.mu(panIdx + 1) - bPan.mu(panIdx)) / deltaL2 +
                           (bPan.mu(panIdx) - bPan.mu(panIdx - 1)) / deltaL1);
            }
            // Transverse
            if (j == 0) { // Forward differencing
                deltaT2 = (bPan.CG.row(panIdx + bPan.nC_) - bPan.CG.row(panIdx)).norm();
                v = (bPan.mu(panIdx + bPan.nC_) - bPan.mu(panIdx)) / deltaT2;
            }
            else if (j == bPan.nS_ - 1) { // Backward differencing
                deltaT1 = (bPan.CG.row(panIdx) - bPan.CG.row(panIdx - bPan.nC_)).norm();
                v = (bPan.mu(panIdx) - bPan.mu(panIdx - bPan.nC_)) / deltaT1;
            }
            else { // Central differencing
                deltaT1 = (bPan.CG.row(panIdx) - bPan.CG.row(panIdx - bPan.nC_)).norm();
                deltaT2 = (bPan.CG.row(panIdx + bPan.nC_) - bPan.CG.row(panIdx)).norm();
                v = 0.5 * ((bPan.mu(panIdx + bPan.nC_) - bPan.mu(panIdx)) / deltaT2 +
                           (bPan.mu(panIdx) - bPan.mu(panIdx - bPan.nC_)) / deltaT1);
            }

            // Transverse to perpendicular velocity
            v = (bPan.t(panIdx, 0) * bPan.p(panIdx, 0)
                 + bPan.t(panIdx, 1) * bPan.p(panIdx, 1)
                 + bPan.t(panIdx, 2) * bPan.p(panIdx, 2)) * v;
            // Normal velocity
            w = bPan.tau(panIdx);

            // Transformation to global axis
            bPan.U(panIdx, 0) = u * bPan.l(panIdx, 0) - v * bPan.p(panIdx, 0) + w * bPan.n(panIdx, 0);
            bPan.U(panIdx, 1) = u * bPan.l(panIdx, 1) - v * bPan.p(panIdx, 1) + w * bPan.n(panIdx, 1);
            bPan.U(panIdx, 2) = u * bPan.l(panIdx, 2) - v * bPan.p(panIdx, 2) + w * bPan.n(panIdx, 2);

            panIdx++;
        }
    }
    // Freestream velocity
    for (int i = 0; i < bPan.nP; ++i) {
        bPan.U(i, 0) += vInf(0);
        bPan.U(i, 1) += vInf(1);
        bPan.U(i, 2) += vInf(2);
    }
    // Source induced velocity
    if (machInf != 0)
        bPan.U += vSigma;
    // Mach number
    if (machInf != 0)
        bPan.M = VectorXd::Zero(bPan.nP);
    else {
        for (int i = 0; i < bPan.nP; ++i)
            bPan.M(i) = bPan.U.row(i).norm()
                        / (1 / (machInf * machInf)
                           + (GAMMA - 1) / 2 - (GAMMA - 1) / 2 * bPan.U.row(i).dot(bPan.U.row(i)));
    }
    //// Aerodynamic forces
    // Pressure coefficient
    if (machInf == 0) {
        for (int i = 0; i < bPan.nP; ++i)
            bPan.cP(i) = 1 - bPan.U.row(i).dot(bPan.U.row(i)) / vInf.dot(vInf);
    }
    else {
        for (int i = 0; i < bPan.nP; ++i)
            bPan.cP(i) = 2 / (GAMMA * machInf * machInf) *
                (pow(1 + (GAMMA - 1) / 2 * machInf * machInf *
                                 (1 - bPan.U.row(i).dot(bPan.U.row(i))), GAMMA / (GAMMA - 1)) - 1);
    }
    // Forces (x,y,z)
    for (int i = 0; i < bPan.nP; ++i) {
        cX += -bPan.cP(i) * bPan.S(i) / sRef * bPan.n(i,0);
        cZ += -bPan.cP(i) * bPan.S(i) / sRef * bPan.n(i,2);
    }
    if (symY) {
        cX = 2*cX;
        cZ = 2*cZ;
    }
    // Forces (l,d,y)
    cD = cos(alpha)*cX + sin(alpha)*cZ;
    cL = -sin(alpha)*cX + cos(alpha)*cZ;

    //// Control display
    cout << "Done!" << endl;
    cout << "Surface velocities: " << bPan.U.rows() << "X" << bPan.U.cols() << endl;
    //for (int i = 0; i < p; ++i)
    //    cout << i << ' ' << vB(i,0) << ' ' << vB(i,1) << ' ' << vB(i,2) << endl;
    cout << "Pressure coefficients: " << bPan.cP.rows() << "X" << bPan.cP.cols() << endl;
    //for (int i = 0; i < p; ++i)
    //    cout << i << ' ' << cP(i) << endl;
    cout << "Aerodynamic forces:" << endl;
    cout << "Cx = " << cX << endl;
    cout << "Cz = " << cZ << endl;
    cout << "Cl = " << cL << endl;
    cout << "Cd = " << cD << endl;
    cout << "L/D = " << cL/cD << endl << endl;
}