//// Compute field variables
// Compute variables (velocity, speed of sound, density, mach number) in the field
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - mgVar: field variabes for minigrid (structure)
// - sp: subpanels indices (structure)
// - dRho : derivative of density
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - mgAIC: body to field AIC for minigrid (structure)
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

void compute_fVars(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Subpanel &sp, MatrixX3d &dRho,
                       Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC, Subpanel_AIC &spAIC) {

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
    // Perturbation velocity (minigrid)
    mgVar.UXbwd.col(0) = mgAIC.BuXbwd * bPan.tau + mgAIC.AuXbwd * bPan.mu + mgAIC.CuXbwd * fPan.sigma;
    mgVar.UXbwd.col(1) = mgAIC.BvXbwd * bPan.tau + mgAIC.AvXbwd * bPan.mu + mgAIC.CvXbwd * fPan.sigma;
    mgVar.UXbwd.col(2) = mgAIC.BwXbwd * bPan.tau + mgAIC.AwXbwd * bPan.mu + mgAIC.CwXbwd * fPan.sigma;
    mgVar.UXfwd.col(0) = mgAIC.BuXfwd * bPan.tau + mgAIC.AuXfwd * bPan.mu + mgAIC.CuXfwd * fPan.sigma;
    mgVar.UXfwd.col(1) = mgAIC.BvXfwd * bPan.tau + mgAIC.AvXfwd * bPan.mu + mgAIC.CvXfwd * fPan.sigma;
    mgVar.UXfwd.col(2) = mgAIC.BwXfwd * bPan.tau + mgAIC.AwXfwd * bPan.mu + mgAIC.CwXfwd * fPan.sigma;
    mgVar.UYbwd.col(0) = mgAIC.BuYbwd * bPan.tau + mgAIC.AuYbwd * bPan.mu + mgAIC.CuYbwd * fPan.sigma;
    mgVar.UYbwd.col(1) = mgAIC.BvYbwd * bPan.tau + mgAIC.AvYbwd * bPan.mu + mgAIC.CvYbwd * fPan.sigma;
    mgVar.UYbwd.col(2) = mgAIC.BwYbwd * bPan.tau + mgAIC.AwYbwd * bPan.mu + mgAIC.CwYbwd * fPan.sigma;
    mgVar.UYfwd.col(0) = mgAIC.BuYfwd * bPan.tau + mgAIC.AuYfwd * bPan.mu + mgAIC.CuYfwd * fPan.sigma;
    mgVar.UYfwd.col(1) = mgAIC.BvYfwd * bPan.tau + mgAIC.AvYfwd * bPan.mu + mgAIC.CvYfwd * fPan.sigma;
    mgVar.UYfwd.col(2) = mgAIC.BwYfwd * bPan.tau + mgAIC.AwYfwd * bPan.mu + mgAIC.CwYfwd * fPan.sigma;
    mgVar.UZbwd.col(0) = mgAIC.BuZbwd * bPan.tau + mgAIC.AuZbwd * bPan.mu + mgAIC.CuZbwd * fPan.sigma;
    mgVar.UZbwd.col(1) = mgAIC.BvZbwd * bPan.tau + mgAIC.AvZbwd * bPan.mu + mgAIC.CvZbwd * fPan.sigma;
    mgVar.UZbwd.col(2) = mgAIC.BwZbwd * bPan.tau + mgAIC.AwZbwd * bPan.mu + mgAIC.CwZbwd * fPan.sigma;
    mgVar.UZfwd.col(0) = mgAIC.BuZfwd * bPan.tau + mgAIC.AuZfwd * bPan.mu + mgAIC.CuZfwd * fPan.sigma;
    mgVar.UZfwd.col(1) = mgAIC.BvZfwd * bPan.tau + mgAIC.AvZfwd * bPan.mu + mgAIC.CvZfwd * fPan.sigma;
    mgVar.UZfwd.col(2) = mgAIC.BwZfwd * bPan.tau + mgAIC.AwZfwd * bPan.mu + mgAIC.CwZfwd * fPan.sigma;

    // Perturbation velocity (with sub-paneling technique)
    for (int jj = 0; jj < sp.sI.size(); jj++) {
    //for (int jj = 0; jj < 1; jj++) {
        int j = sp.sI[jj]; // panel global index
        int n = j % bPan.nC_; // panel chordwise index
        int m = j / bPan.nC_; // panel spanwise index
        // Interpolation of singularities from centers to vertices
        vSing = interp_ctv(j, n, m, bPan);
        // Interpolation of singularities from vertices to sub-panel
        sInterp = interp_sp(j, bPan, sp, vSing(0,0), vSing(1,0), vSing(2,0), vSing(3,0),
                            vSing(0,1), vSing(0,1), vSing(2,1), vSing(3,1));

        for (int ii = 0; ii < sp.fI[jj].size(); ii++) {
            int i = sp.fI[jj][ii]; // cell index
            // Velocity induces by each sub-panel
            for (int k = 0; k < sp.NS; k++) {
                // Cell center
                fPan.U(i,0) += spAIC.Au[jj](ii,k) * sInterp(k,0) + spAIC.Bu[jj](ii,k) * sInterp(k,1);
                fPan.U(i,1) += spAIC.Av[jj](ii,k) * sInterp(k,0) + spAIC.Bv[jj](ii,k) * sInterp(k,1);
                fPan.U(i,2) += spAIC.Aw[jj](ii,k) * sInterp(k,0) + spAIC.Bw[jj](ii,k) * sInterp(k,1);
                // Minigrid
                mgVar.UXbwd(i,0) += spAIC.AuXbwd[jj](ii,k) * sInterp(k,0) + spAIC.BuXbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UXbwd(i,1) += spAIC.AvXbwd[jj](ii,k) * sInterp(k,0) + spAIC.BvXbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UXbwd(i,2) += spAIC.AwXbwd[jj](ii,k) * sInterp(k,0) + spAIC.BwXbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UXfwd(i,0) += spAIC.AuXfwd[jj](ii,k) * sInterp(k,0) + spAIC.BuXfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UXfwd(i,1) += spAIC.AvXfwd[jj](ii,k) * sInterp(k,0) + spAIC.BvXfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UXfwd(i,2) += spAIC.AwXfwd[jj](ii,k) * sInterp(k,0) + spAIC.BwXfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYbwd(i,0) += spAIC.AuYbwd[jj](ii,k) * sInterp(k,0) + spAIC.BuYbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYbwd(i,1) += spAIC.AvYbwd[jj](ii,k) * sInterp(k,0) + spAIC.BvYbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYbwd(i,2) += spAIC.AwYbwd[jj](ii,k) * sInterp(k,0) + spAIC.BwYbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYfwd(i,0) += spAIC.AuYfwd[jj](ii,k) * sInterp(k,0) + spAIC.BuYfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYfwd(i,1) += spAIC.AvYfwd[jj](ii,k) * sInterp(k,0) + spAIC.BvYfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UYfwd(i,2) += spAIC.AwYfwd[jj](ii,k) * sInterp(k,0) + spAIC.BwYfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZbwd(i,0) += spAIC.AuZbwd[jj](ii,k) * sInterp(k,0) + spAIC.BuZbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZbwd(i,1) += spAIC.AvZbwd[jj](ii,k) * sInterp(k,0) + spAIC.BvZbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZbwd(i,2) += spAIC.AwZbwd[jj](ii,k) * sInterp(k,0) + spAIC.BwZbwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZfwd(i,0) += spAIC.AuZfwd[jj](ii,k) * sInterp(k,0) + spAIC.BuZfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZfwd(i,1) += spAIC.AvZfwd[jj](ii,k) * sInterp(k,0) + spAIC.BvZfwd[jj](ii,k) * sInterp(k,1);
                mgVar.UZfwd(i,2) += spAIC.AwZfwd[jj](ii,k) * sInterp(k,0) + spAIC.BwZfwd[jj](ii,k) * sInterp(k,1);
            }
        }
    }

    // Freestream component
    for (int i = 0; i < fPan.nE; ++i) {
        fPan.U.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UXbwd.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UXfwd.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UYbwd.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UYfwd.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UZbwd.row(fPan.eIdx(i)) += vInf.transpose();
        mgVar.UZfwd.row(fPan.eIdx(i)) += vInf.transpose();
    }

    // TODO: check if useful
    // Wake treatment
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        if (fPan.wMap(idx)) {
            fPan.U(idx, 0) = 0.5 * (fPan.U(idx - fPan.nX, 0) + fPan.U(idx + fPan.nX, 0));
            mgVar.UXbwd(idx, 0) = 0.5 * (mgVar.UXbwd(idx - fPan.nX, 0) + mgVar.UXbwd(idx + fPan.nX, 0));
            mgVar.UXfwd(idx, 0) = 0.5 * (mgVar.UXfwd(idx - fPan.nX, 0) + mgVar.UXfwd(idx + fPan.nX, 0));
            mgVar.UYbwd(idx, 0) = 0.5 * (mgVar.UYbwd(idx - fPan.nX, 0) + mgVar.UYbwd(idx + fPan.nX, 0));
            mgVar.UYfwd(idx, 0) = 0.5 * (mgVar.UYfwd(idx - fPan.nX, 0) + mgVar.UYfwd(idx + fPan.nX, 0));
            mgVar.UZbwd(idx, 0) = 0.5 * (mgVar.UZbwd(idx - fPan.nX, 0) + mgVar.UZbwd(idx + fPan.nX, 0));
            mgVar.UZfwd(idx, 0) = 0.5 * (mgVar.UZfwd(idx - fPan.nX, 0) + mgVar.UZfwd(idx + fPan.nX, 0));

            fPan.U(idx, 1) = 0.5 * (fPan.U(idx - fPan.nX, 1) + fPan.U(idx + fPan.nX, 1));
            mgVar.UXbwd(idx, 1) = 0.5 * (mgVar.UXbwd(idx - fPan.nX, 1) + mgVar.UXbwd(idx + fPan.nX, 1));
            mgVar.UXfwd(idx, 1) = 0.5 * (mgVar.UXfwd(idx - fPan.nX, 1) + mgVar.UXfwd(idx + fPan.nX, 1));
            mgVar.UYbwd(idx, 1) = 0.5 * (mgVar.UYbwd(idx - fPan.nX, 1) + mgVar.UYbwd(idx + fPan.nX, 1));
            mgVar.UYfwd(idx, 1) = 0.5 * (mgVar.UYfwd(idx - fPan.nX, 1) + mgVar.UYfwd(idx + fPan.nX, 1));
            mgVar.UZbwd(idx, 1) = 0.5 * (mgVar.UZbwd(idx - fPan.nX, 1) + mgVar.UZbwd(idx + fPan.nX, 1));
            mgVar.UZfwd(idx, 1) = 0.5 * (mgVar.UZfwd(idx - fPan.nX, 1) + mgVar.UZfwd(idx + fPan.nX, 1));

            fPan.U(idx, 2) = 0;
            mgVar.UXbwd(idx, 2) = 0;
            mgVar.UXfwd(idx, 2) = 0;
            mgVar.UYbwd(idx, 2) = 0;
            mgVar.UYfwd(idx, 2) = 0;
            mgVar.UZbwd(idx, 2) = 0;
            mgVar.UZfwd(idx, 2) = 0;
        } else
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
    // Density for minigrid
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        mgVar.rhoXbwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UXbwd.row(idx).dot(mgVar.UXbwd.row(idx))),
                1 / (GAMMA - 1));
        mgVar.rhoXfwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UXfwd.row(idx).dot(mgVar.UXfwd.row(idx))),
                1 / (GAMMA - 1));
        mgVar.rhoYbwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UYbwd.row(idx).dot(mgVar.UYbwd.row(idx))),
                1 / (GAMMA - 1));
        mgVar.rhoYfwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UYfwd.row(idx).dot(mgVar.UYfwd.row(idx))),
                1 / (GAMMA - 1));
        mgVar.rhoZbwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UZbwd.row(idx).dot(mgVar.UZbwd.row(idx))),
                1 / (GAMMA - 1));
        mgVar.rhoZfwd(idx) = pow(
                1 + (GAMMA - 1) / 2 * Minf * Minf * (1 - mgVar.UZfwd.row(idx).dot(mgVar.UZfwd.row(idx))),
                1 / (GAMMA - 1));
    }
}