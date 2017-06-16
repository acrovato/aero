//// Solve field
// Compute field variables and field sources
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
#include "solve_field.h"
#include "compute_fVars.h"

using namespace std;
using  namespace Eigen;

#define GAMMA 1.4
#define CMU 1.0
#define M_C 1.0

void solve_field(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Subpanel &sp,
                 Body_AIC &b2fAIC, Field2field_AIC &f2fAIC, Subpanel_AIC &spAIC) {

    int idx; // counter
    int ost = 0; // cell offset for finite differencing
    double mu, deltaSigmaX, deltaSigmaY, deltaSigmaZ; // Variables for artificial viscosity
    double dU, dV, dW; // Variables for residual

    //// Begin
    cout << "Computing field variables... " << flush;

    //// Field variables
    compute_fVars(Minf, vInf, bPan, fPan, sp, b2fAIC, f2fAIC, spAIC);

    //// Field sources
    // Density gradient
    idx = 0;
    for (int j = 0; j < fPan.nY; ++j) {
        for (int k = 0; k < fPan.nZ; ++k) {
            for (int i = 0; i < fPan.nX; ++i) {

                if (!fPan.fMap(idx)) {
                    fPan.dRho(idx, 0) = 0;
                    fPan.dRho(idx, 1) = 0;
                    fPan.dRho(idx, 2) = 0;
                    dU = 0;
                    dV = 0;
                    dW = 0;
                    fPan.sigma(idx) = 0;
                    fPan.epsilon(idx) = 0;
                }
                else {
                    // X-derivative
                    ost = 1;
                    if (i == 0) {
                        if (fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 0) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaX;
                            dU = (fPan.U(idx + ost, 0) - fPan.U(idx, 0)) / fPan.deltaX;}
                        else {
                            fPan.dRho(idx, 0) = 0;
                            dU = 0; }
                    } else if (i == fPan.nX - 1) {
                        if (fPan.fMap(idx - ost)) {
                            fPan.dRho(idx, 0) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaX;
                            dU = (fPan.U(idx, 0) - fPan.U(idx - ost, 0)) / fPan.deltaX;}
                        else {
                            fPan.dRho(idx, 0) = 0;
                            dU = 0;}
                    } else {
                        if (fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 0) = 0.5 * (fPan.rho(idx + ost) - fPan.rho(idx - ost)) / fPan.deltaX;
                            dU = 0.5 * (fPan.U(idx + ost,0) - fPan.U(idx - ost,0)) / fPan.deltaX;}
                        else if (fPan.fMap(idx - ost) && !fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 0) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaX;
                            dU = (fPan.U(idx,0) - fPan.U(idx - ost,0)) / fPan.deltaX;}
                        else if (!fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 0) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaX;
                            dU = (fPan.U(idx + ost,0) - fPan.U(idx,0)) / fPan.deltaX;}
                        else {
                            fPan.dRho(idx, 0) = 0;
                            dU = 0;}
                    }
                    // Y-derivative
                    ost = fPan.nX*fPan.nZ;
                    if (j == 0) {
                        if (fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 1) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaY;
                            dV = (fPan.U(idx + ost,1) - fPan.U(idx,1)) / fPan.deltaY;}
                        else {
                            fPan.dRho(idx, 1) = 0;
                            dV = 0;}
                    } else if (j == fPan.nY - 1) {
                        if (fPan.fMap(idx - ost)) {
                            fPan.dRho(idx, 1) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaY;
                            dV = (fPan.U(idx,1) - fPan.U(idx - ost,1)) / fPan.deltaY;}
                        else {
                            fPan.dRho(idx, 1) = 0;
                            dV = 0;}
                    } else {
                        if (fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 1) = 0.5 * (fPan.rho(idx + ost) - fPan.rho(idx - ost)) / fPan.deltaY;
                            dV = 0.5 * (fPan.U(idx + ost,1) - fPan.U(idx - ost,1)) / fPan.deltaY;}
                        else if (fPan.fMap(idx - ost) && !fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 1) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaY;
                            dV = (fPan.U(idx,1) - fPan.U(idx - ost,1)) / fPan.deltaY;}
                        else if (!fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 1) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaY;
                            dV = (fPan.U(idx + ost,1) - fPan.U(idx,1)) / fPan.deltaY;}
                        else {
                            fPan.dRho(idx, 1) = 0;
                            dV = 0;}
                    }
                    // Z-derivative
                    ost = fPan.nX;
                    if (k == 0) {
                        if (fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 2) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaZ;
                            dW = (fPan.U(idx + ost,2) - fPan.U(idx,2)) / fPan.deltaZ;}
                        else {
                            fPan.dRho(idx, 2) = 0;
                            dW = 0;}
                    } else if (k == fPan.nZ - 1) {
                        if (fPan.fMap(idx - ost)) {
                            fPan.dRho(idx, 2) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaZ;
                            dW = (fPan.U(idx,2) - fPan.U(idx - ost,2)) / fPan.deltaZ;}
                        else {
                            fPan.dRho(idx, 2) = 0;
                            dW = 0;}
                    } else {
                        if (fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 2) = 0.5 * (fPan.rho(idx + ost) - fPan.rho(idx - ost)) / fPan.deltaZ;
                            dW = 0.5 * (fPan.U(idx + ost,2) - fPan.U(idx - ost,2)) / fPan.deltaZ;}
                        else if (fPan.fMap(idx - ost) && !fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 2) = (fPan.rho(idx) - fPan.rho(idx - ost)) / fPan.deltaZ;
                            dW = (fPan.U(idx,2) - fPan.U(idx - ost,2)) / fPan.deltaZ;}
                        else if (!fPan.fMap(idx - ost) && fPan.fMap(idx + ost)) {
                            fPan.dRho(idx, 2) = (fPan.rho(idx + ost) - fPan.rho(idx)) / fPan.deltaZ;
                            dW = (fPan.U(idx + ost,2) - fPan.U(idx,2)) / fPan.deltaZ;}
                        else {
                            fPan.dRho(idx, 2) = 0;
                            dW = 0;}
                    }
                    // Field source
                    fPan.sigma(idx) = -1 / (fPan.rho(idx)) * fPan.U.row(idx).dot(fPan.dRho.row(idx));
                    // Residual
                    fPan.epsilon(idx) = dU + dV + dW - fPan.sigma(idx);
                }
                idx++;
            }
        }
    }

    // TODO 1) Artificial density (pros: physical; cons: does not work on MG or RG/MG)
    // TODO 2) Artificial viscosity (pros: works on RG/MG; cons: cut through surface, not physical)
    // TODO NB) With current form, x-upwinding gives same results as s-upwinding.

    // TODO Implement safeguard for when cell is on border on computational domain (especially on y)
    //// Artificial viscosity
    idx = 0;
    for (int j = 0; j < fPan.nY; ++j) {
        for (int k = 0; k < fPan.nZ; ++k) {
            for (int i = 0; i < fPan.nX; ++i) {

                if (fPan.fMap(idx) && fPan.M(idx) > M_C) {
                    mu = CMU * (1 - M_C * M_C / (fPan.M(idx) * fPan.M(idx)));

                    // X-contribution
                    ost = 1;
                    if (fPan.U(idx, 0) > 0 && fPan.fMap(idx - ost) && i != 0)
                        deltaSigmaX = fPan.sigma(idx) - fPan.sigma(idx - ost);
                    else if (fPan.U(idx, 0) < 0 && fPan.fMap(idx + ost) && i != fPan.nX-1)
                        deltaSigmaX = fPan.sigma(idx + ost) - fPan.sigma(idx);
                    else
                        deltaSigmaX = 0;
                    // Y-contribution
                    ost = fPan.nX*fPan.nZ;
                    if (fPan.U(idx, 1) > 0 && fPan.fMap(idx - ost) && j != 0)
                        deltaSigmaY = fPan.sigma(idx) - fPan.sigma(idx - ost);
                    else if (fPan.U(idx, 1) < 0 && fPan.fMap(idx + ost) && j != fPan.nY-1)
                        deltaSigmaY = fPan.sigma(idx + ost) - fPan.sigma(idx);
                    else
                        deltaSigmaY = 0;
                    // Z-contribution
                    ost = fPan.nX;
                    if (fPan.U(idx, 2) > 0 && fPan.fMap(idx - ost) && k != 0)
                        deltaSigmaZ = fPan.sigma(idx) - fPan.sigma(idx - ost);
                    else if (fPan.U(idx, 2) < 0 && fPan.fMap(idx + ost) && k != fPan.nZ-1)
                        deltaSigmaZ = fPan.sigma(idx + ost) - fPan.sigma(idx);
                    else
                        deltaSigmaZ = 0;

                    fPan.sigma(idx) -= mu / fPan.U.row(idx).norm() *
                                       (fPan.U(idx, 0) * deltaSigmaX + fPan.U(idx, 1) * deltaSigmaY + fPan.U(idx, 2) * deltaSigmaZ);
                }
                idx ++;
            }
        }
    }

    cout << "Done!" << endl;
    cout << "Max. Mach number: " << fPan.M.maxCoeff() << endl;
    cout << "Field sources strength: " << fPan.sigma.norm() << endl << endl;
}