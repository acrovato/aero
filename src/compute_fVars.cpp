//// Compute field variables
// Compute variables (velocity, speed of sound, density, mach number) in the field
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - mgVar: field variabes for minigrid (structure)
// - dRho : derivative of density
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - mgAIC: body to field AIC for minigrid (structure)

#include <iostream>
#include <Eigen/Dense>
#include "compute_fVars.h"

#define GAMMA 1.4

using namespace std;
using namespace Eigen;

void compute_fVars(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, MatrixX3d &dRho,
                       Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC) {

    int idx; // temporary counter

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