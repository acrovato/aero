//// Solver
// Build AIC matrices and solve iteratively the "body" and the "field"
// Body module: compute impermeability condition and solve linear system of equations (corresponding to the
// Linear Potential Equation)
// Field module: compute through AIC and FD the variables in the field, as well as the field sources (source term in
// the Full Potential Equation)
//
// I/O:
// - numC: list of numerical parameters (structure)
// - symY: defines symmetry about Y axis
// - sRef: reference surface of the full wing
// - alpha: freestream angle of attack
// - vInf: freestream velocity vector
// - Minf: freestream Mach number
// - bPan: (network of) body panels (structure)
// - wPan: (network of) wake panels (structure)
// - fPan: field panels (structure)
// - sp: sub-panels (structure)
// - cL: lift coefficient
// - cD: drag coefficient
//
// Output:
// Outputs:
// - 0, if function succeeded
// - 1, if function failed

#include <iostream>
#include <Eigen/Dense>
#include "solver.h"
#include "id_subpanel.h"
#include "build_AIC.h"
#include "solve_body.h"
#include "solve_field.h"
#include "compute_sVars.h"

#define ANSI_COLOR_RED     "\x1b[1;31m"
#define ANSI_COLOR_GREEN   "\x1b[1;32m"
#define ANSI_COLOR_YELLOW  "\x1b[1;33m"
#define ANSI_COLOR_BLUE    "\x1b[1;34m"
#define ANSI_COLOR_MAGENTA "\x1b[1;35m"
#define ANSI_COLOR_CYAN    "\x1b[1;36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define NDIM 3

using namespace std;
using namespace Eigen;

int solver(Numerical_CST &numC, bool symY, double sRef, double alpha, Vector3d &vInf, double Minf,
           Network &bPan, Network &wPan, Field &fPan, Subpanel &sp, double &cL, double &cD) {


    //// Begin solver
    cout << "*********************" << endl;
    cout << "*Beginning solver...*" << endl;
    cout << "*********************" << endl;
    cout << endl;

    //// Initialization
    // AIC matrices
    Body_AIC b2bAIC = {}, b2fAIC = {}; // body to body and body to field
    Field2field_AIC f2fAIC = {}; // field to field
    Field2body_AIC f2bAIC = {}; // field to field and field to body
    Subpanel_AIC spAIC = {}; // body to field (sub-panel)

    // Singularities
    bPan.tau.resize(bPan.nP);
    bPan.mu.resize(bPan.nP);
    fPan.sigma = VectorXd::Zero(fPan.nF);
    // Flow variables // TODO check if initialization is useful
    fPan.phi = VectorXd::Zero(fPan.nF);
    fPan.M = VectorXd::Zero(fPan.nF);
    fPan.U = MatrixX3d::Zero(fPan.nF, NDIM);
    fPan.rho = VectorXd::Zero(fPan.nF);
    fPan.dRho = MatrixX3d::Zero(fPan.nF, NDIM);
    fPan.a = VectorXd::Zero(fPan.nF);
    // Numerics
    fPan.epsilon = VectorXd::Zero(fPan.nF);


    // Temporary variables
    int itCnt = 0; // Global iteration counter
    double deltaSigma0 = 0; // Initial sigma change
    VectorXd sigmaTmp, deltaSigma; // To store sigma and delta sigma during iteration
    deltaSigma.resize(fPan.nF);
    sigmaTmp.resize(fPan.nF);
    VectorXd RHS; // Right hand side
    RHS.resize(bPan.nP);

    MatrixX3d vSigma; // Velocity induced by field sources on body
    vSigma = MatrixXd::Zero(bPan.nP, NDIM);

    //// Identify sub-panels
    id_subpanel(bPan, fPan, sp, spAIC);

    //// Build AIC matrices
    build_AIC(symY, bPan, wPan, fPan, b2bAIC, b2fAIC, f2fAIC, f2bAIC, sp, spAIC);

    //// Solver
    // Field Panel Method
    if (Minf != 0) {
        cout << "|¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯|" << endl;
        cout << "|Compressible computation.|" << endl;
        cout << " ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ " << endl << endl;
        do {
            // Store new field sources
            sigmaTmp = fPan.sigma;
            // Panel prediction and boundary condition
            solve_body(vInf, RHS, vSigma, bPan, b2bAIC);
            // Field correction
            solve_field(Minf, vInf, bPan, fPan, sp, b2fAIC, f2fAIC, spAIC);
            // Source induced velocity (to recompute B.C.)
            vSigma.col(0) = f2bAIC.Cu * fPan.sigma;
            vSigma.col(1) = f2bAIC.Cv * fPan.sigma;
            vSigma.col(2) = f2bAIC.Cw * fPan.sigma;
            // Stop criterion
            for (int i = 0; i < fPan.nF; ++i)
                deltaSigma(i) = abs(fPan.sigma(i) - sigmaTmp(i));
            if (!itCnt)
                deltaSigma0 = deltaSigma.norm();
            if (isnan(deltaSigma.norm())) {
                cout << ANSI_COLOR_RED << "∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨" << endl;
                cout << ">> Process diverged at iteration #" << itCnt + 1 << "!" << endl;
                cout << "∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧" << ANSI_COLOR_RESET << endl << endl;
                exit(EXIT_FAILURE);
            }
            cout << "Relative source change at iteration " << itCnt << ": " << log10(deltaSigma.norm()/deltaSigma0) << endl;
            cout << "FPE global residual at iteration " << itCnt << ": " << log10(fPan.epsilon.norm()) << endl << endl;
            itCnt++;
        } while(log10(deltaSigma.norm()/deltaSigma0) > -numC.RRED);
        //TODO Consider using true residual (div(U) - sigma -> 0). Currently impossible since accuracy.
        cout << ANSI_COLOR_GREEN << "∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨" << endl;
        cout << ">> Process converged in " << itCnt << " iteration(s)!" << endl;
        cout << "∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧" << ANSI_COLOR_RESET << endl << endl;
    }
    // Panel Method
    else {
        cout << "|¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯|" << endl;
        cout << "|Incompressible computation.|" << endl;
        cout << " ¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯ " << endl << endl;
        solve_body(vInf, RHS, vSigma, bPan, b2bAIC);
    }

    // Surface velocity and pressure computation
    compute_sVars(symY, sRef, alpha, Minf, vInf, vSigma, bPan, cL, cD);

    //// End solver
    cout << "********************" << endl;
    cout << "*Solver successful!*" << endl;
    cout << "********************" << endl;
    return 0;
}