//// Solver
// Build AIC matrices and solve iteratively the "body" and the "field"
// Body module: compute impermeability condition and solve linear system of equations (corresponding to the
// Linear Potential Equation)
// Field module: compute through AIC and FD the variables in the field, as well as the field sources (source term in
// the Full Potential Equation)
//
// I/O:
// - symY: defines symmetry about Y axis
// - sRef: reference surface of the full wing
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - wPan: (network of) wake panels (structure)
// - fPan: field panels (structure)
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
    Body_AIC b2bAIC = {}; // body to body
    Body2field_AIC b2fAIC = {}; // body to field
    Field_AIC f2fAIC, f2bAIC = {}; // field to field and field to body
    Minigrid_AIC mgAIC = {}; // field to field and field to body (minigrid)
    Subpanel_AIC spAIC = {}; // body to field (sub-panel)

    // Singularities
    bPan.tau.resize(bPan.nP);
    bPan.mu.resize(bPan.nP);
    fPan.sigma = VectorXd::Zero(fPan.nF);
    fPan.sigmaTilda = VectorXd::Zero(fPan.nF);
    fPan.dSigma = MatrixX3d::Zero(fPan.nF, NDIM);
    // Flow variables // TODO check if initialization is useful
    fPan.M = VectorXd::Zero(fPan.nF);
    fPan.U = MatrixX3d::Zero(fPan.nF, NDIM);
    fPan.rho = VectorXd::Zero(fPan.nF);
    fPan.dRho = MatrixX3d::Zero(fPan.nF, NDIM);
    fPan.a = VectorXd::Zero(fPan.nF);

    // Temporary variables
    int itCnt = 0; // Global iteration counter
    double resInit = 0; // Initial residual
    VectorXd sigmaTmp, res; // To store sigma during iteration and residual
    res.resize(fPan.nF);
    sigmaTmp.resize(fPan.nF);
    VectorXd RHS; // Right hand side
    RHS.resize(bPan.nP);

    MatrixX3d vSigma; // Velocity induced by field sources on body
    vSigma = MatrixXd::Zero(bPan.nP, NDIM);
    Minigrid mgVar = {}; // Minigrid variables
    mgVar.rhoXbwd.resize(fPan.nF);
    mgVar.rhoXfwd.resize(fPan.nF);
    mgVar.rhoYbwd.resize(fPan.nF);
    mgVar.rhoYfwd.resize(fPan.nF);
    mgVar.rhoZbwd.resize(fPan.nF);
    mgVar.rhoZfwd.resize(fPan.nF);
    mgVar.UXbwd.resize(fPan.nF, NDIM);
    mgVar.UXfwd.resize(fPan.nF, NDIM);
    mgVar.UYbwd.resize(fPan.nF, NDIM);
    mgVar.UYfwd.resize(fPan.nF, NDIM);
    mgVar.UZbwd.resize(fPan.nF, NDIM);
    mgVar.UZfwd.resize(fPan.nF, NDIM);
    mgVar.rhoXmidbwd.resize(fPan.nF);
    mgVar.rhoXmidfwd.resize(fPan.nF);
    mgVar.rhoYmidbwd.resize(fPan.nF);
    mgVar.rhoYmidfwd.resize(fPan.nF);
    mgVar.rhoZmidbwd.resize(fPan.nF);
    mgVar.rhoZmidfwd.resize(fPan.nF);
    mgVar.UXmidbwd.resize(fPan.nF, NDIM);
    mgVar.UXmidfwd.resize(fPan.nF, NDIM);
    mgVar.UYmidbwd.resize(fPan.nF, NDIM);
    mgVar.UYmidfwd.resize(fPan.nF, NDIM);
    mgVar.UZmidbwd.resize(fPan.nF, NDIM);
    mgVar.UZmidfwd.resize(fPan.nF, NDIM);
    mgVar.dRhoXmidbwd.resize(fPan.nF, NDIM);
    mgVar.dRhoXmidfwd.resize(fPan.nF, NDIM);
    mgVar.dRhoYmidbwd.resize(fPan.nF, NDIM);
    mgVar.dRhoYmidfwd.resize(fPan.nF, NDIM);
    mgVar.dRhoZmidbwd.resize(fPan.nF, NDIM);
    mgVar.dRhoZmidfwd.resize(fPan.nF, NDIM);
    mgVar.sigmaXmidbwd.resize(fPan.nF);
    mgVar.sigmaXmidfwd.resize(fPan.nF);
    mgVar.sigmaYmidbwd.resize(fPan.nF);
    mgVar.sigmaYmidfwd.resize(fPan.nF);
    mgVar.sigmaZmidbwd.resize(fPan.nF);
    mgVar.sigmaZmidfwd.resize(fPan.nF);

    //// Identify sub-panels
    id_subpanel(bPan, fPan, sp, spAIC);

    //// Build AIC matrices
    build_AIC(symY, bPan, wPan, fPan, b2bAIC, b2fAIC, f2fAIC, f2bAIC, mgAIC, sp, spAIC);

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
            solve_field(Minf, vInf, bPan, fPan, mgVar, sp, b2fAIC, f2fAIC, mgAIC, spAIC);
            // Source induced velocity (to recompute B.C.)
            vSigma.col(0) = f2bAIC.Cu * fPan.sigma;
            vSigma.col(1) = f2bAIC.Cv * fPan.sigma;
            vSigma.col(2) = f2bAIC.Cw * fPan.sigma;
            // Stop criterion
            for (int i = 0; i < fPan.nF; ++i)
                res(i) = abs(fPan.sigma(i) - sigmaTmp(i));
            if (!itCnt)
                resInit = res.norm();
            if (isnan(res.norm())) {
                cout << "∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨" << endl;
                cout << ">>Process diverged at iteration #" << itCnt + 1 << "!<<" << endl;
                cout << "∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧" << endl << endl;
                exit(EXIT_FAILURE);
            }
            cout << "Rel. residual at iteration " << itCnt << ": " << log10(res.norm()/resInit) << endl;
            cout << "Max. residual at iteration " << itCnt << ": " << log10(res.maxCoeff()) << endl << endl;
            itCnt++;
        } while(itCnt < 15);// (log10(res.norm()/resInit) > -numC.RRED);
        cout << "∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨∨" << endl;
        cout << ">>Process converged in " << itCnt << " iteration(s)!<<" << endl;
        cout << "∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧∧" << endl << endl;
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