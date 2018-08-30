//// Main
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
// - sp: sub-panels (structure)
//
// *Physics*
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - vInf: freestream velocity vector
// - cL: lift coefficient
// - cD: drag coefficient
//
// *Numerics*
// - numC: numerical parameters (structure)

#include <iostream>
#include <Eigen/Dense>
#include <ctime>
#include "pre.h"
#include "solver.h"
#include "post.h"

#define ANSI_COLOR_BLUE    "\x1b[1;34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;
using namespace Eigen;

int main( int argc, char *argv[] ) {

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
    // Hello World
    time_t now = time(0); // get time now
    char* localNow = ctime(&now);
    cout << ANSI_COLOR_BLUE;
    cout << "***********************************" << endl;
    cout << "**              \\_/              **" << endl;
    cout << "**     \\_______O(_)O_______/     **" << endl;
    cout << "**         _                     **" << endl;
    cout << "**        / \\   __  __  __       **" << endl;
    cout << "**       / _ \\ |__||__||  |      **" << endl;
    cout << "**      /_/ \\_\\|__ |  \\|__|      **" << endl;
    cout << "***********************************";
    cout << ANSI_COLOR_RESET << endl;
    cout << "Hi! My name is Aero v1.0-1809" << endl;
    cout << "Solver started on " << localNow << endl;

    // Check parameters
    if (argc != 3) {
        cout << "Incorrect number of parameters provided!" << endl;
        cout << "Usage: ./aero <pathToConfigFile(.cfg)> <pathToGridgridFile(.pts)>." << endl;
        return 1;
    }

    // Set time counter
    clock_t start = clock();
    clock_t startS = clock();

    // Pre-processing
    pre(argv, numC, symY, sRef, Minf, alpha, vInf, bPan, wPan, fPan, sp);
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