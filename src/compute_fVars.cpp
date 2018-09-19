//// Field variables computation
// Compute variables (velocity, speed of sound, density, mach number) in the field
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - sp: sub-panel (structure)
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - spAIC: body to field AIC for subpanels (structure)

/* Copyright (C) 2018 Adrien Crovato */

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
    for (int i = 0; i < fPan.nE; i++) {
        int f = fPan.eIdx(i);

        // X-derivative
        double fb=0, ff=0;
        if (fPan.fbdMap(f,0) && fPan.fwdMap(f,0))
            fb = (fPan.phi(f) - fPan.phi(f - 1)) / fPan.deltaX;
        if (fPan.fbdMap(f,1) && fPan.fwdMap(f,1))
            ff = (fPan.phi(f + 1) - fPan.phi(f)) / fPan.deltaX;
        if (fPan.fbdMap(f,0) && fPan.fwdMap(f,0) && fPan.fbdMap(f,1) && fPan.fwdMap(f,1))
            fPan.U(f,0) = 0.5*(fb+ff);
        else
            fPan.U(f,0) = fb+ff;

        // Y-derivative
        fb=0; ff=0;
        if (fPan.fbdMap(f,2) && fPan.fwdMap(f,2))
            fb = (fPan.phi(f) - fPan.phi(f - fPan.nX*fPan.nZ)) / fPan.deltaY;
        if (fPan.fbdMap(f,3) && fPan.fwdMap(f,3))
            ff = (fPan.phi(f + fPan.nX*fPan.nZ) - fPan.phi(f)) / fPan.deltaY;
        if (fPan.fbdMap(f,2) && fPan.fwdMap(f,2) && fPan.fbdMap(f,3) && fPan.fwdMap(f,3))
            fPan.U(f,1) = 0.5*(fb+ff);
        else
            fPan.U(f,1) = fb+ff;

        // Z-derivative
        fb=0; ff=0;
        if (fPan.fbdMap(f,4) && fPan.fwdMap(f,4))
            fb = (fPan.phi(f) - fPan.phi(f - fPan.nX)) / fPan.deltaZ;
        if (fPan.fbdMap(f,5) && fPan.fwdMap(f,5))
            ff = (fPan.phi(f + fPan.nX) - fPan.phi(f)) / fPan.deltaZ;
        if (fPan.fbdMap(f,4) && fPan.fwdMap(f,4) && fPan.fbdMap(f,5) && fPan.fwdMap(f,5))
            fPan.U(f,2) = 0.5*(fb+ff);
        else
            fPan.U(f,2) = fb+ff;
    }
    // Freestream component
    for (int i = 0; i < fPan.nE; ++i) {
        fPan.U.row(fPan.eIdx(i)) += vInf.transpose();
    }

    //// Thermodynamic variables
    for (int i = 0; i < fPan.nE; ++i) {
        int f = fPan.eIdx(i);
        // Speed of sound
        fPan.a(f) = sqrt(
                1/(Minf * Minf) + (GAMMA - 1) / 2 - (GAMMA - 1) / 2 * fPan.U.row(f).dot(fPan.U.row(f)));
        // Mach number
        fPan.M(f) = fPan.U.row(f).norm() / fPan.a(f);
        // Density
        fPan.rho(f) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - fPan.U.row(f).dot(fPan.U.row(f))),
                1 / (GAMMA - 1));
    }
}