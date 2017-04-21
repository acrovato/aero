//// Solve body
// Compute sources on the body surface from the updated  impermeability B.C. and solve for doublets on the body surface
//
// Ref: Katz & Plotkin (2001), Low-Speed Aerodynamics
// 1) tau = - (U_inf + Usigma) \cdot n
// 2) mu = A \ (B * tau)
//
// I/O:
// - vInf: freestream velocity vector
// - RHS: right-hand side of the linear system of equations
// - vSigma: field source induced velocty
// - bPan: (network of) body panels (structure)
// - b2bAIC: body to body AIC (structure)

#include <iostream>
#include <Eigen/Dense>
#include "solve_body.h"

using namespace std;
using namespace Eigen;

void solve_body(Vector3d &vInf, VectorXd &RHS, MatrixX3d &vSigma, Network &bPan, Body_AIC &b2bAIC) {

    //// Begin
    cout << "Computing surface singularities... " << flush;

    //// Update B.C.
    // BC: tau_i = n_i * V_inf + n_i * V_sigma
    for (int i = 0; i < bPan.nP; ++i)
        bPan.tau(i) = - bPan.n.row(i).dot(vInf.transpose()) - bPan.n.row(i).dot(vSigma.row(i));

    //// Solve
    // Compute RHS
    RHS = - b2bAIC.B * bPan.tau;
    // Solve: A*mu + B*sigma = 0 (Dirichlet: interior potential = 0)
    #ifdef ON_UNIX
        bPan.mu = b2bAIC.A.householderQr().solve(RHS);
    #else
        panDoublet = muAIC.fullPivLu().solve(RHS);
    #endif

    //// Control display
    cout << "Done!" << endl;
    cout << "Surface sources strength: " << bPan.tau.norm() << endl;
    //cout << "Sources min. value: " << panSource.minCoeff() << endl;
    //cout << "Sources max. value: " << panSource.maxCoeff() << endl;
    //for (int i = 0; i < nP; ++i)
    //    cout << i << ' ' << panSource(i) << endl;
    cout << "Surface doublets strength: " << bPan.mu.norm() << endl;
    //cout << "Doublets min. value: " << panDoublet.minCoeff() << endl;
    //cout << "Doublets max. value: " << panDoublet.maxCoeff() << endl;
    //for (int i = 0; i < nP; ++i)
    //    cout << i << ' ' << panDoublet(i) << endl;
    cout << "System solved with relative error: " << (b2bAIC.A * bPan.mu - RHS).norm() / RHS.norm() << endl;
    cout << endl;
}
