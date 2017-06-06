//// main file
// Link program modules (pre, solver, post)
//
// Variables accessed throughout the code
//
// *Geometry*
// - symY: defines symmetry about Y axis
// - sRef: reference surface of the full wing
//
// *Networks and Field*
// - bPan: (network of) body panels (structure)
// - wPan: (network of) wake panels (structure)
// - fPan: field panels (structure)
//
// *Physics*
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - vInf: freestream velocity vector
// - cL: lift coefficient
// - cD: drag coefficient

//// Acronym
// - STREAM: Swift Transonic Enhanced Aerodynamic Modeling
//		JET-STREAM, STEAM
// - QtAero: Quick Aerodynamics
// - STAMPD: Swift Transonic Aerodynamic Modeling for Preliminary Design
// - STAR(PAD): Swift Transonic Aerodynamics (for Preliminary Aircraft Design)
// - SAMPAD: Swift Aerodynamic Modeling for Preliminary Aircraft Design
// - CELIA: C - Enhanced Lightweight Aerodynamics
// - CeLiA: C Library for Aerodynamics
// - CeLIA: C Library for swIft (or quIck) Aerodynamics
// - CLELiA: C Library for Enhanced Linear Aerodynamics
// - ST(R)ELA: Swift Transonic Enhanced Lightweight Aerodynamics
// - SCALP: Swift Computational Aerodynamic Loads Prediction
// - AeroMAD: Aerodynamic Modeling for Aircraft Design
// - ETA: Enhanced Transonic Aerodynamics

#include <iostream>
#include <Eigen/Dense>
#include "pre.h"
#include "solver.h"
#include "post.h"

using namespace std;
using namespace Eigen;

int main() {

    //// Variable definition
	// Geometry
    bool symY = 0;
    double sRef = 1;
    // Constants
    Numerical_CST numC = {};
	// Freestream
    double Minf = 0;
    double alpha = 0;
    Vector3d vInf(1.0, 0.0, 0.0);
	// Surface and wake panels
    // TODO: if several Networks are considered, use std::vector <network>
	Network bPan = {};
	Network wPan = {};
	// Field cells
	Field fPan = {};
    // Sub-panels
    Subpanel sp = {};
	// Forces
    double cL = 0, cD = 0;

    //// Begin flow
    // Hello World: Fast and Accurate Steady Transonic Aerodynamics for AeroElastic Tailoring
    time_t now = time(0); // get time now
    char* localNow = ctime(&now);
    cout << "**********************************" << endl;
    cout << "**             \\_/              **" << endl;
    cout << "**     \\______O(_)O_______/     **" << endl;
    cout << "**                              **" << endl;
    cout << "**   ____     _     _    __     **" << endl;
    cout << "**  /  __|__ | |   (_)  /  \\    **" << endl;
    cout << "** (  |__|__|| |__ | | / /\\ \\   **" << endl;
    cout << "**  \\____|__ |____||_|/_/¯¯\\_\\  **" << endl;
    cout << "**********************************" << endl << endl;
    cout << "Hi! My name is CeLiA v1.0-1704" << endl;
    cout << "Solver started on " << localNow << endl;

    // Set time counter
    clock_t start = clock();
    clock_t startS = clock();

    // Pre-processing
    pre(numC, symY, sRef, Minf, alpha, vInf, bPan, wPan, fPan, sp);
    clock_t endS = clock();
    cout << "Preprocessing time: " << (endS - startS) / (double) CLOCKS_PER_SEC << "s" << endl << endl;

    // Solver
    startS = clock();
    solver(numC, symY, sRef, alpha, vInf, Minf, bPan, wPan, fPan, sp, cL, cD);
    endS = clock();
    cout << "Solver time: " << (endS - startS) / (double) CLOCKS_PER_SEC << "s" << endl << endl;

    // Post-processing
    startS = clock();
    post(sRef, alpha, Minf, bPan, fPan, cL, cD);
    endS = clock();
    cout << "Postprocessing time: " << (endS - startS) / (double) CLOCKS_PER_SEC << "s" << endl << endl;

    // Total time
    clock_t end = clock();
    cout << "***** Run summary *****" << endl;
    cout << "Total time: " << (end - start) / (double) CLOCKS_PER_SEC << "s" << endl << endl;

    return 0;
}