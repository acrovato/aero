//// Compressible module (field solver)
// Compute field variables and field sources
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - sp: sub-panels (structure)
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - spAIC: body to field AIC for subpanels (structure)

/* Copyright (C) 2018 Adrien Crovato */

#include <iostream>
#include <Eigen/Dense>
#include "solve_field.h"
#include "compute_fVars.h"

using namespace std;
using  namespace Eigen;

#define GAMMA 1.4
#define CMU 1.0
#define M_C 1.0

void solve_field(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Subpanel &sp,
                 Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Subpanel_AIC &spAIC) {

    //// Begin
    cout << "Computing field variables... " << flush;

    //// Field variables
    compute_fVars(Minf, vInf, bPan, fPan, sp, b2fAIC, f2fAIC, spAIC);

    //// Field sources
    // Density gradient
    for (int i = 0; i < fPan.nE; i++) {
        int f = fPan.eIdx(i);
        double dU, dV, dW; // Variables for residual

        // X-derivative
        double fb=0, ff=0, fb2=0, ff2=0;
        if (fPan.fbdMap(f,0)) {
            fb = (fPan.rho(f) - fPan.rho(f - 1)) / fPan.deltaX;
            fb2 = (fPan.U(f,0) - fPan.U(f - 1,0)) / fPan.deltaX;
        }
        if (fPan.fbdMap(f,1)) {
            ff = (fPan.rho(f + 1) - fPan.rho(f)) / fPan.deltaX;
            ff2 = (fPan.U(f + 1,0) - fPan.U(f,0)) / fPan.deltaX;
        }
        if (fPan.fbdMap(f,0) && fPan.fbdMap(f,1)) {
            fPan.dRho(f,0) = 0.5*(fb+ff);
            dU = 0.5*(fb2+ff2);
        }
        else {
            fPan.dRho(f,0) = fb+ff;
            dU = fb2+ff2;
        }

        // Y-derivative
        fb=0, ff=0, fb2=0, ff2=0;
        if (fPan.fbdMap(f,2)) {
            fb = (fPan.rho(f) - fPan.rho(f - fPan.nX*fPan.nZ)) / fPan.deltaY;
            fb2 = (fPan.U(f,1) - fPan.U(f - fPan.nX*fPan.nZ,1)) / fPan.deltaY;
        }
        if (fPan.fbdMap(f,3)) {
            ff = (fPan.rho(f + fPan.nX*fPan.nZ) - fPan.rho(f)) / fPan.deltaY;
            ff2 = (fPan.U(f + fPan.nX*fPan.nZ,1) - fPan.U(f,1)) / fPan.deltaY;
        }
        if (fPan.fbdMap(f,2) && fPan.fbdMap(f,3)) {
            fPan.dRho(f,1) = 0.5*(fb+ff);
            dV = 0.5*(fb2+ff2);
        }
        else {
            fPan.dRho(f,1) = fb+ff;
            dV = fb2+ff2;
        }

        // Z-derivative
        fb=0, ff=0, fb2=0, ff2=0;
        if (fPan.fbdMap(f,4)) {
            fb = (fPan.rho(f) - fPan.rho(f - fPan.nX)) / fPan.deltaZ;
            fb2 = (fPan.U(f,2) - fPan.U(f - fPan.nX,2)) / fPan.deltaZ;
        }
        if (fPan.fbdMap(f,5)) {
            ff = (fPan.rho(f + fPan.nX) - fPan.rho(f)) / fPan.deltaZ;
            ff2 = (fPan.U(f + fPan.nX,2) - fPan.U(f,2)) / fPan.deltaZ;
        }
        if (fPan.fbdMap(f,4) && fPan.fbdMap(f,5)) {
            fPan.dRho(f,2) = 0.5*(fb+ff);
            dW = 0.5*(fb2+ff2);
        }
        else {
            fPan.dRho(f,2) = fb+ff;
            dW = fb2+ff2;
        }

        // Field source
        fPan.sigma(f) = -1 / (fPan.rho(f)) * fPan.U.row(f).dot(fPan.dRho.row(f));
        // Residual
        fPan.epsilon(f) = dU + dV + dW - fPan.sigma(f);
    }

    // TODO 1) Artificial density (pros: physical; cons: does not work on MG or RG/MG)
    // TODO 2) Artificial viscosity (pros: works on RG/MG; cons: cut through surface, not physical)
    // TODO NB) With current form, x-upwinding gives same results as s-upwinding.

    //// Artificial viscosity
    for (int i = 0; i < fPan.nE; i++) {
        int f = fPan.eIdx(i);
        double mu, deltaSigmaX, deltaSigmaY, deltaSigmaZ; // Variables for artificial viscosity

        if (fPan.M(f) > M_C) {
            mu = CMU * (1 - M_C * M_C / (fPan.M(f) * fPan.M(f)));

            // X-contribution
            if (fPan.U(f,0) > 0 && fPan.fbdMap(f,0))
                deltaSigmaX = fPan.sigma(f) - fPan.sigma(f - 1);
            else if (fPan.U(f,0) < 0 && fPan.fbdMap(f,1))
                deltaSigmaX = fPan.sigma(f + 1) - fPan.sigma(f);
            else
                deltaSigmaX = 0;
            // Y-contribution
            if (fPan.U(f,1) > 0 && fPan.fbdMap(f,2))
                deltaSigmaY = fPan.sigma(f) - fPan.sigma(f - fPan.nX*fPan.nZ);
            else if (fPan.U(f,1) < 0 && fPan.fbdMap(f,3))
                deltaSigmaY = fPan.sigma(f + fPan.nX*fPan.nZ) - fPan.sigma(f);
            else
                deltaSigmaY = 0;
            // Z-contribution
            if (fPan.U(f,2) > 0 && fPan.fbdMap(f,4))
                deltaSigmaZ = fPan.sigma(f) - fPan.sigma(f - fPan.nX);
            else if (fPan.U(f,2) < 0 && fPan.fbdMap(f,5))
                deltaSigmaZ = fPan.sigma(f + fPan.nX) - fPan.sigma(f);
            else
                deltaSigmaZ = 0;

            fPan.sigma(f) -= mu / fPan.U.row(f).norm() *
                               (fPan.U(f, 0) * deltaSigmaX + fPan.U(f, 1) * deltaSigmaY + fPan.U(f, 2) * deltaSigmaZ);
        }
    }

    cout << "Done!" << endl;
    cout << "Max. Mach number: " << fPan.M.maxCoeff() << endl;
    cout << "Field sources strength: " << fPan.sigma.norm() << endl << endl;
}