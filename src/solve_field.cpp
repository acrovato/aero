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
#define M_C 1.0

void solve_field(double Minf, Vector3d &vInf, Network &bPan, Field &fPan, Minigrid &mgVar, Subpanel &sp,
                 Body2field_AIC &b2fAIC, Field_AIC &f2fAIC, Minigrid_AIC &mgAIC, Subpanel_AIC &spAIC) {

    int idx = 0; // counter

    //// Begin
    cout << "Computing field variables... " << flush;

    //// Field variables
    compute_fVars(Minf, vInf, bPan, fPan, mgVar, sp, b2fAIC, f2fAIC, mgAIC, spAIC);

    //// Field sources
    // Velocity and density interpolation
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Velociy
        mgVar.UXmidbwd.row(idx) = 0.5 * (mgVar.UXbwd.row(idx) + fPan.U.row(idx));
        mgVar.UXmidfwd.row(idx) = 0.5 * (fPan.U.row(idx) + mgVar.UXfwd.row(idx));
        mgVar.UYmidbwd.row(idx) = 0.5 * (mgVar.UYbwd.row(idx) + fPan.U.row(idx));
        mgVar.UYmidfwd.row(idx) = 0.5 * (fPan.U.row(idx) + mgVar.UYfwd.row(idx));
        mgVar.UZmidbwd.row(idx) = 0.5 * (mgVar.UZbwd.row(idx) + fPan.U.row(idx));
        mgVar.UZmidfwd.row(idx) = 0.5 * (fPan.U.row(idx) + mgVar.UZfwd.row(idx));
        // Density
        mgVar.rhoXmidbwd(idx) = 0.5 * (mgVar.rhoXbwd(idx) + fPan.rho(idx));
        mgVar.rhoXmidfwd(idx) = 0.5 * (fPan.rho(idx) + mgVar.rhoXfwd(idx));
        mgVar.rhoYmidbwd(idx) = 0.5 * (mgVar.rhoYbwd(idx) + fPan.rho(idx));
        mgVar.rhoYmidfwd(idx) = 0.5 * (fPan.rho(idx) + mgVar.rhoYfwd(idx));
        mgVar.rhoZmidbwd(idx) = 0.5 * (mgVar.rhoZbwd(idx) + fPan.rho(idx));
        mgVar.rhoZmidfwd(idx) = 0.5 * (fPan.rho(idx) + mgVar.rhoZfwd(idx));
    }

    // Density gradient
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Central point
        // Subsonic point, central differencing
        //if (fPan.M(idx) < M_C) {
            fPan.dRho(idx, 0) = 0.5 * (mgVar.rhoXfwd(idx) - mgVar.rhoXbwd(idx)) / fPan.deltaMG;
            fPan.dRho(idx, 1) = 0.5 * (mgVar.rhoYfwd(idx) - mgVar.rhoYbwd(idx)) / fPan.deltaMG;
            fPan.dRho(idx, 2) = 0.5 * (mgVar.rhoZfwd(idx) - mgVar.rhoZbwd(idx)) / fPan.deltaMG;
        //}
        // Supersonic point, upwind differencing
        //else {
        //    if (fPan.U(idx, 0) > 0)
        //        fPan.dRho(idx, 0) = (fPan.rho(idx) - mgVar.rhoXbwd(idx)) / fPan.deltaMG;
        //    else
        //        fPan.dRho(idx, 0) = (mgVar.rhoXfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
        //    if (fPan.U(idx, 1) > 0)
        //        fPan.dRho(idx, 1) = (fPan.rho(idx) - mgVar.rhoYbwd(idx)) / fPan.deltaMG;
        //    else
        //        fPan.dRho(idx, 1) = (mgVar.rhoYfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
        //    if (fPan.U(idx, 2) > 0)
        //        fPan.dRho(idx, 2) = (fPan.rho(idx) - mgVar.rhoZbwd(idx)) / fPan.deltaMG;
        //    else
        //        fPan.dRho(idx, 2) = (mgVar.rhoZfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
        //}
        // Minigrid midpoint
        mgVar.dRhoXmidbwd(idx,0) = (fPan.rho(idx) - mgVar.rhoXbwd(idx)) / fPan.deltaMG;
        mgVar.dRhoXmidbwd(idx,1) = 0.5*((mgVar.rhoXbwd(idx)+mgVar.rhoYfwd(idx)) - (mgVar.rhoXbwd(idx)+mgVar.rhoYbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoXmidbwd(idx,2) = 0.5*((mgVar.rhoXbwd(idx)+mgVar.rhoZfwd(idx)) - (mgVar.rhoXbwd(idx)+mgVar.rhoZbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoXmidfwd(idx,0) = (mgVar.rhoXfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
        mgVar.dRhoXmidfwd(idx,1) = 0.5*((mgVar.rhoXfwd(idx)+mgVar.rhoYfwd(idx)) - (mgVar.rhoXfwd(idx)+mgVar.rhoYbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoXmidfwd(idx,2) = 0.5*((mgVar.rhoXfwd(idx)+mgVar.rhoZfwd(idx)) - (mgVar.rhoXfwd(idx)+mgVar.rhoZbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoYmidbwd(idx,0) = 0.5*((mgVar.rhoYbwd(idx)+mgVar.rhoXfwd(idx)) - (mgVar.rhoYbwd(idx)+mgVar.rhoXbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoYmidbwd(idx,1) = (fPan.rho(idx) - mgVar.rhoYbwd(idx)) / fPan.deltaMG;
        mgVar.dRhoYmidbwd(idx,2) = 0.5*((mgVar.rhoYbwd(idx)+mgVar.rhoZfwd(idx)) - (mgVar.rhoYbwd(idx)+mgVar.rhoZbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoYmidfwd(idx,0) = 0.5*((mgVar.rhoYfwd(idx)+mgVar.rhoXfwd(idx)) - (mgVar.rhoYfwd(idx)+mgVar.rhoXbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoYmidfwd(idx,1) = (mgVar.rhoYfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
        mgVar.dRhoYmidfwd(idx,2) = 0.5*((mgVar.rhoYfwd(idx)+mgVar.rhoZfwd(idx)) - (mgVar.rhoYfwd(idx)+mgVar.rhoZbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoZmidbwd(idx,0) = 0.5*((mgVar.rhoZbwd(idx)+mgVar.rhoXfwd(idx)) - (mgVar.rhoZbwd(idx)+mgVar.rhoXbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoZmidbwd(idx,1) = 0.5*((mgVar.rhoZbwd(idx)+mgVar.rhoYfwd(idx)) - (mgVar.rhoZbwd(idx)+mgVar.rhoYbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoZmidbwd(idx,2) = (fPan.rho(idx) - mgVar.rhoZbwd(idx)) / fPan.deltaMG;
        mgVar.dRhoZmidfwd(idx,0) = 0.5*((mgVar.rhoZfwd(idx)+mgVar.rhoXfwd(idx)) - (mgVar.rhoZfwd(idx)+mgVar.rhoXbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoZmidfwd(idx,1) = 0.5*((mgVar.rhoZfwd(idx)+mgVar.rhoYfwd(idx)) - (mgVar.rhoZfwd(idx)+mgVar.rhoYbwd(idx))) / fPan.deltaMG;
        mgVar.dRhoZmidfwd(idx,2) = (mgVar.rhoZfwd(idx) - fPan.rho(idx)) / fPan.deltaMG;
    }

    // Field sources
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Central point
        fPan.sigma(idx) = -1 / fPan.rho(idx) * fPan.U.row(idx).dot(fPan.dRho.row(idx));
        // Minigrid midpoint
        mgVar.sigmaXmidbwd(idx) = -1 / mgVar.rhoXmidbwd(idx) * mgVar.UXmidbwd.row(idx).dot(mgVar.dRhoXmidbwd.row(idx));
        mgVar.sigmaXmidfwd(idx) = -1 / mgVar.rhoXmidfwd(idx) * mgVar.UXmidfwd.row(idx).dot(mgVar.dRhoXmidfwd.row(idx));
        mgVar.sigmaYmidbwd(idx) = -1 / mgVar.rhoYmidbwd(idx) * mgVar.UYmidbwd.row(idx).dot(mgVar.dRhoYmidbwd.row(idx));
        mgVar.sigmaYmidfwd(idx) = -1 / mgVar.rhoYmidfwd(idx) * mgVar.UYmidfwd.row(idx).dot(mgVar.dRhoYmidfwd.row(idx));
        mgVar.sigmaZmidbwd(idx) = -1 / mgVar.rhoZmidbwd(idx) * mgVar.UZmidbwd.row(idx).dot(mgVar.dRhoZmidbwd.row(idx));
        mgVar.sigmaZmidfwd(idx) = -1 / mgVar.rhoZmidfwd(idx) * mgVar.UZmidfwd.row(idx).dot(mgVar.dRhoZmidfwd.row(idx));
    }

    // Artificial viscosity
    for (int i = 0; i < fPan.nE; ++i) {
        idx = fPan.eIdx(i);
        // Field sources gradient
        fPan.dSigma(idx,0) = (mgVar.sigmaXmidfwd(idx) - mgVar.sigmaXmidbwd(idx)) / fPan.deltaMG;
        fPan.dSigma(idx,1) = (mgVar.sigmaYmidfwd(idx) - mgVar.sigmaYmidbwd(idx)) / fPan.deltaMG;
        fPan.dSigma(idx,2) = (mgVar.sigmaZmidfwd(idx) - mgVar.sigmaZmidbwd(idx)) / fPan.deltaMG;
        // Viscosity
        fPan.sigmaTilda(idx) = (1 / (fPan.M(idx)*fPan.M(idx)) - 1) / fPan.U.row(idx).norm()
                               * (fPan.U(idx,0) * fPan.dSigma(idx,0) * fPan.deltaMG/2
                                  + fPan.U(idx,1) * fPan.dSigma(idx,1) * fPan.deltaMG/2
                                  + fPan.U(idx,2) * fPan.dSigma(idx,2) * fPan.deltaMG/2);
        // Field source update
        if (fPan.M(idx) > 1)
            fPan.sigma(idx) += fPan.sigmaTilda(idx);
    }

    cout << "Done!" << endl;
    cout << "Max. Mach number: " << fPan.M.maxCoeff() << endl;
    cout << "Field sources strength: " << fPan.sigma.norm() << endl << endl;
}