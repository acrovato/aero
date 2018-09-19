//// Pre-processing
// Load settings and geometry from external files *.cfg and *.pts (located in /IO/)
// Store information into relevant structures
//
// Inputs:
// - arvg: command line argument provided to main (contains path to config and grid files)
// - numC: list of numerical parameters (structure)
// - symY: defines symmetry about Y axis
// - sRef: reference surface of the full wing
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - vInf: freestream velocity vector
// - bPan: (network of) body panels (structure)
// - wPan: (network of) wake panels (structure)
// - fPan: field panels (structure)
// Outputs:
// - 0, if function succeeded
// - 1, if function failed

/* Copyright (C) 2018 Adrien Crovato */

#include <iostream>
#include <array>
#include <Eigen/Dense>

#include "pre.h"
#include "read_config.h"
#include "read_sgrid.h"
#include "create_panel.h"
#include "create_wake.h"
#include "create_field.h"
#include "map_field.h"
#include "map_derivatives.h"

#define ANSI_COLOR_BLUE    "\x1b[1;34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;
using namespace Eigen;

int pre(char *argv[], Numerical_CST &numC, bool &symY, double &sRef, double &machInf, double &AoA, Vector3d &vInf,
        Network &bPan, Network &wPan, Field &fPan, Subpanel &sp) {

    //// Initialization
    // Temporary arrays to store grid information
    array<array<double, 3>, 8> box{}; // Corner points defining domain
    MatrixX3d sGrid ; // Surface grid corner points

    // Path definition
    string configParamPath = argv[1];
    string surfGridPath = argv[2];

    //// Begin preprocessing
    cout << ANSI_COLOR_BLUE;
    cout << "*****************************" << endl;
    cout << "*Beginning pre-processing...*" << endl;
    cout << "*****************************";
    cout << ANSI_COLOR_RESET << endl;

    // Read config and grid from files
    read_config(configParamPath, numC, symY, sRef, machInf, AoA, box, fPan, sp);
    read_sgrid(surfGridPath, sGrid, bPan);

    // Set freestream velocity and speed of sound
    vInf(0) = cos(AoA);
    vInf(1) = 0;
    vInf(2) = sin(AoA);

    // Create surface grid
    create_panel(sGrid, bPan);
    create_wake(sGrid, bPan, wPan);

    // Create volume grid, map cells and derivatives
    create_field(box, fPan);
    map_field(sGrid, numC, bPan, fPan);
    map_derivatives(sGrid, numC, bPan, wPan, fPan);

    //// End preprocessing
    cout << ANSI_COLOR_BLUE;
    cout << "****************************" << endl;
    cout << "*Pre-processing successful!*" << endl;
    cout << "****************************";
    cout << ANSI_COLOR_RESET << endl;
    return 0;
}