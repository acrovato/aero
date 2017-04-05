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
// - dRho : derivative of density
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
#define M_C 1.0

void solve_field(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Subpanel &sp, MatrixX3d &dRho,
                 Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC, Subpanel_AIC &spAIC) {

    int idx; // counter

    //// Begin
    cout << "Computing field variables... " << flush;

    //// Field variables
    compute_fVars(Minf, vInf, bPan, fPan, mgVar, sp, dRho, b2fAIC, f2fAIC, mgAIC, spAIC);

    //// Field sources
    // X-derivative of density
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Subsonic point, central differencing
        dRho(idx, 0) = 0.5 * (mgVar.rhoXfwd(idx) - mgVar.rhoXbwd(idx)) / fPan.deltaMG;
    }
    // Y-derivative of density
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Subsonic point, central differencing
        dRho(idx, 1) = 0.5 * (mgVar.rhoYfwd(idx) - mgVar.rhoYbwd(idx)) / fPan.deltaMG;
    }
    // Z-derivative of density
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Subsonic point, central differencing
        dRho(idx, 2) = 0.5 * (mgVar.rhoZfwd(idx) - mgVar.rhoZbwd(idx)) / fPan.deltaMG;
    }

    // Field sources
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        fPan.sigma(idx) = -1 / fPan.rho(idx) * fPan.U.row(idx).dot(dRho.row(idx));
    }

    cout << "Done!" << endl;
    cout << "Max. Mach number: " << fPan.M.maxCoeff() << endl;
    cout << "Field sources strength: " << fPan.sigma.norm() << endl << endl;
}