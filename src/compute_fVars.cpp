//// Compute field variables
// Compute variables (velocity, speed of sound, density, mach number) in the field
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - sp: subpanels indices (structure)
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - spAIC: body to field AIC for subpanels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "compute_fVars.h"
#include "interp_ctv.h"
#include "interp_sp.h"

#define GAMMA 1.4
#define NV 4 // number of vertices

using namespace std;
using namespace Eigen;

void compute_fVars(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Subpanel &sp,
                       Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Subpanel_AIC &spAIC) {

    //// Initialization
    // Temporary variables
    int idx; // counter
    MatrixXd vSing; // corner interpolated singularities
    vSing.resize(NV,2);
    MatrixXd sInterp; // sub-panel interoplated singularities
    sInterp.resize(sp.NS,2);

    //// Velocity
    // Perturbation velocity
    fPan.U.col(0) = b2fAIC.Bu * bPan.tau + b2fAIC.Au * bPan.mu + f2fAIC.Cu * fPan.sigma;
    fPan.U.col(1) = b2fAIC.Bv * bPan.tau + b2fAIC.Av * bPan.mu + f2fAIC.Cv * fPan.sigma;
    fPan.U.col(2) = b2fAIC.Bw * bPan.tau + b2fAIC.Aw * bPan.mu + f2fAIC.Cw * fPan.sigma;

    // Perturbation velocity (with sub-paneling technique)
    for (int jj = 0; jj < sp.sI.size(); jj++) {
    //for (int jj = 0; jj < 1; jj++) {
        int j = sp.sI[jj]; // panel global index
        int n = j % bPan.nC_; // panel chordwise index
        int m = j / bPan.nC_; // panel spanwise index
        // Interpolation of singularities from centers to vertices
        vSing = interp_ctv(j, n, m, bPan);
        // Interpolation of singularities from vertices to sub-panel
        sInterp = interp_sp(j, bPan, sp,
                            vSing(0,0), vSing(1,0), vSing(2,0), vSing(3,0),
                            vSing(0,1), vSing(1,1), vSing(2,1), vSing(3,1));

        for (int ii = 0; ii < sp.fI[jj].size(); ii++) {
            int i = sp.fI[jj][ii]; // cell index
            // Velocity induces by each sub-panel
            for (int k = 0; k < sp.NS; k++) {
                // Cell center
                fPan.U(i,0) += spAIC.Au[jj](ii,k) * sInterp(k,0) + spAIC.Bu[jj](ii,k) * sInterp(k,1);
                fPan.U(i,1) += spAIC.Av[jj](ii,k) * sInterp(k,0) + spAIC.Bv[jj](ii,k) * sInterp(k,1);
                fPan.U(i,2) += spAIC.Aw[jj](ii,k) * sInterp(k,0) + spAIC.Bw[jj](ii,k) * sInterp(k,1);
            }
        }
    }

    // Freestream component
    for (int i = 0; i < fPan.nE; ++i) {
        fPan.U.row(fPan.eIdx(i)) += vInf.transpose();
    }

    // TODO Needs improvement
    // TODO a) Use same algorithm in wake and field to avoid cutting though surface
    // TODO b) use cutoff distance between cell and to wake to use (or not) extrapolation between k+1 and k+2 for cell k
    // Wake treatment
    for (int i = 0; i < fPan.nW; ++i) {
        idx = fPan.wIdx(i);
        if (!fPan.wMap(idx-fPan.nX) && !fPan.wMap(idx+fPan.nX)) {
            fPan.U(idx, 0) = 0.5 * (fPan.U(idx - fPan.nX, 0) + fPan.U(idx + fPan.nX, 0));
            fPan.U(idx, 1) = 0.5 * (fPan.U(idx - fPan.nX, 1) + fPan.U(idx + fPan.nX, 1));
            fPan.U(idx, 2) = 0;
        }
        else if (!fPan.wMap(idx-fPan.nX) && fPan.wMap(idx+fPan.nX)) {
            fPan.U(idx, 0) = fPan.U(idx - fPan.nX, 0);
            fPan.U(idx, 1) = fPan.U(idx - fPan.nX, 1);
            fPan.U(idx, 2) = 0;
        }
        else if (fPan.wMap(idx-fPan.nX) && !fPan.wMap(idx+fPan.nX)) {
            fPan.U(idx, 0) = fPan.U(idx + fPan.nX, 0);
            fPan.U(idx, 1) = fPan.U(idx + fPan.nX, 1);
            fPan.U(idx, 2) = 0;
        }
        else
            continue;
    }

    //// Thermodynamic variables
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Speed of sound
        fPan.a(idx) = sqrt(
                1/(Minf * Minf) + (GAMMA - 1) / 2 - (GAMMA - 1) / 2 * fPan.U.row(idx).dot(fPan.U.row(idx)));
        // Mach number
        fPan.M(idx) = fPan.U.row(idx).norm() / fPan.a(idx);
        // Density
        fPan.rho(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - fPan.U.row(idx).dot(fPan.U.row(idx))),
                1 / (GAMMA - 1));
    }
}