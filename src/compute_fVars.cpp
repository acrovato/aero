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
                       Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Subpanel_AIC &spAIC) {

    //// Initialization
    // Temporary variables
    int idx; // counter
    int ost = 0; // cell offset for finite differencing
    MatrixXd vSing; // corner interpolated singularities
    vSing.resize(NV,2);
    MatrixXd sInterp; // sub-panel interoplated singularities
    sInterp.resize(sp.NS,2);

    //// Potential
    // Perturbation potential
    fPan.phi = b2fAIC.B * bPan.tau + b2fAIC.A * bPan.mu + f2fAIC.C * fPan.sigma;

    // Perturbation potential (with sub-paneling technique)
    // TODO Subpaneling induces small assymmetry. Check why?
    for (int jj = 0; jj < sp.sI.size(); jj++) {
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
                fPan.phi(i) += spAIC.A[jj](ii,k) * sInterp(k,0) + spAIC.B[jj](ii,k) * sInterp(k,1);
            }
        }
    }

    //// Velocity
    // Perturbation velocity
    idx = 0;
    for (int j = 0; j < fPan.nY; ++j) {
        for (int k = 0; k < fPan.nZ; ++k) {
            for (int i = 0; i < fPan.nX; ++i) {

                if (!fPan.fMap(idx)) {
                    fPan.U(idx, 0) = 0;
                    fPan.U(idx, 1) = 0;
                    fPan.U(idx, 2) = 0;
                }
                else {
                    // X-derivative
                    ost = 1;
                    if (i == 0) {
                        if (fPan.fMap(idx + ost))
                            fPan.U(idx, 0) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaX;
                        else
                            fPan.U(idx, 0) = 0;
                    } else if (i == fPan.nX - 1) {
                        if (fPan.fMap(idx - ost))
                            fPan.U(idx, 0) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaX;
                        else
                            fPan.U(idx, 0) = 0;
                    } else {
                        if (fPan.fMap(idx - ost) && fPan.fMap(idx + ost))
                            fPan.U(idx, 0) = 0.5 * (fPan.phi(idx + ost) - fPan.phi(idx - ost)) / fPan.deltaX;
                        else if (fPan.fMap(idx - ost) && !fPan.fMap(idx + ost))
                            fPan.U(idx, 0) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaX;
                        else if (!fPan.fMap(idx - ost) && fPan.fMap(idx + ost))
                            fPan.U(idx, 0) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaX;
                        else
                            fPan.U(idx, 0) = 0;
                    }
                    // Y-derivative
                    ost = fPan.nX*fPan.nZ;
                    if (j == 0) {
                        if (fPan.fMap(idx + ost))
                            fPan.U(idx, 1) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaY;
                        else
                            fPan.U(idx, 1) = 0;
                    } else if (j == fPan.nY - 1) {
                        if (fPan.fMap(idx - ost))
                            fPan.U(idx, 1) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaY;
                        else
                            fPan.U(idx, 1) = 0;
                    } else {
                        if (fPan.fMap(idx - ost) && fPan.fMap(idx + ost))
                            fPan.U(idx, 1) = 0.5 * (fPan.phi(idx + ost) - fPan.phi(idx - ost)) / fPan.deltaY;
                        else if (fPan.fMap(idx - ost) && !fPan.fMap(idx + ost))
                            fPan.U(idx, 1) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaY;
                        else if (!fPan.fMap(idx - ost) && fPan.fMap(idx + ost))
                            fPan.U(idx, 1) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaY;
                        else
                            fPan.U(idx, 1) = 0;
                    }
                    // Z-derivative
                    // TODO Add derivative through wake treatment, similar to derivative through surface treatment
                    if (fPan.wMap(idx))
                        fPan.U(idx, 2) = 0;
                    else {
                        ost = fPan.nX;
                        if (k == 0) {
                            if (fPan.fMap(idx + ost) && !fPan.wMap(idx + ost))
                                fPan.U(idx, 2) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaZ;
                            else
                                fPan.U(idx, 2) = 0;
                        } else if (k == fPan.nZ - 1) {
                            if (fPan.fMap(idx - ost) && !fPan.wMap(idx - ost))
                                fPan.U(idx, 2) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaZ;
                            else
                                fPan.U(idx, 2) = 0;
                        } else {
                            if ((fPan.fMap(idx - ost) && !fPan.wMap(idx - ost)) &&
                                (fPan.fMap(idx + ost) && !fPan.wMap(idx + ost)))
                                fPan.U(idx, 2) = 0.5 * (fPan.phi(idx + ost) - fPan.phi(idx - ost)) / fPan.deltaZ;
                            else if ((fPan.fMap(idx - ost) && !fPan.wMap(idx - ost)) &&
                                     (!fPan.fMap(idx + ost) || fPan.wMap(idx + ost)))
                                fPan.U(idx, 2) = (fPan.phi(idx) - fPan.phi(idx - ost)) / fPan.deltaZ;
                            else if ((!fPan.fMap(idx - ost) || fPan.wMap(idx - ost)) &&
                                     (fPan.fMap(idx + ost) && !fPan.wMap(idx + ost)))
                                fPan.U(idx, 2) = (fPan.phi(idx + ost) - fPan.phi(idx)) / fPan.deltaZ;
                            else
                                fPan.U(idx, 2) = 0;
                        }
                    }
                }
                idx++;
            }
        }
    }
    // Freestream component
    for (int i = 0; i < fPan.nE; ++i) {
        fPan.U.row(fPan.eIdx(i)) += vInf.transpose();
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