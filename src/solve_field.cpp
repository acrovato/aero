//// Solve field
// Compute field variables and field sources
//
// I/O:
// - Minf: freestream Mach number
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - fPan: field panels (structure)
// - mgVar: field variabes for minigrid (structure)
// - sp: subpanels indices (structure)
// - b2fAIC: body to field AIC (structure)
// - f2fAIC: field to field AIC (structure)
// - mgAIC: body to field AIC for minigrid (structure)
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

void solve_field(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Subpanel &sp,
                 Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC, Subpanel_AIC &spAIC) {

    int idx = 0; // counter
    double mu, deltaSigmaX, deltaSigmaY, deltaSigmaZ;

    //// Begin
    cout << "Computing field variables... " << flush;

    //// Field variables
    compute_fVars(Minf, vInf, bPan, fPan, mgVar, sp, b2fAIC, f2fAIC, mgAIC, spAIC);

    //// Field sources
    // Density gradient
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Subsonic point, central differencing
        fPan.dRho(idx, 0) = 0.5 * (mgVar.rhoXfwd(idx) - mgVar.rhoXbwd(idx)) / fPan.deltaMG;
        fPan.dRho(idx, 1) = 0.5 * (mgVar.rhoYfwd(idx) - mgVar.rhoYbwd(idx)) / fPan.deltaMG;
        fPan.dRho(idx, 2) = 0.5 * (mgVar.rhoZfwd(idx) - mgVar.rhoZbwd(idx)) / fPan.deltaMG;

        // Field sources
        fPan.sigma(idx) = -1 / (fPan.rho(idx)) * fPan.U.row(idx).dot(fPan.dRho.row(idx));
    }

    // TODO 1) Artificial density (pros: physical; cons: does not work on MG or RG/MG)
    // TODO 2) Artificial viscosity (pros: works on RG/MG; cons: cut through surface, not physical)
    // TODO NB) With current form, x-upwinding gives same results as s-upwinding.

    // TODO Implement safeguard (cells not all defined on RG!!!)
    //// Artificial viscosity
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);

        if (fPan.M(idx) > M_C) {
            mu = CMU * (1 - M_C * M_C / (fPan.M(idx) * fPan.M(idx)));
        if (fPan.U(idx, 0) > 0)
            deltaSigmaX = fPan.sigma(idx) - fPan.sigma(idx - 1);
        else
            deltaSigmaX = fPan.sigma(idx + 1) - fPan.sigma(idx);
        if (fPan.U(idx, 1) > 0 && idx > fPan.nX * fPan.nZ)
            deltaSigmaY = fPan.sigma(idx) - fPan.sigma(idx - fPan.nX * fPan.nZ);
        else if (fPan.U(idx, 1) < 0 && fPan.nX * fPan.nZ * (fPan.nY - 1))
            deltaSigmaY = fPan.sigma(idx + fPan.nX * fPan.nZ) - fPan.sigma(idx);
        else
            deltaSigmaY = 0;
        if (fPan.U(idx, 2) > 0)
            deltaSigmaZ = fPan.sigma(idx) - fPan.sigma(idx - fPan.nX);
        else
            deltaSigmaZ = fPan.sigma(idx + fPan.nX) - fPan.sigma(idx);

        fPan.sigma(idx) -= mu / fPan.U.row(idx).norm() *
                           (fPan.U(idx, 0) * deltaSigmaX + fPan.U(idx, 1) * deltaSigmaY + fPan.U(idx, 2) * deltaSigmaZ);
        }
    }

    cout << "Done!" << endl;
    cout << "Max. Mach number: " << fPan.M.maxCoeff() << endl;
    cout << "Field sources strength: " << fPan.sigma.norm() << endl << endl;
}