//// Preprocessing
// Load settings and geometry from external files *.cfg and *.pts (located in /IO/)
// Store information into relevant structures
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
// Outputs:
// - 0, if function succeeded
// - 1, if function failed

#include <iostream>
#include <array>
#include <Eigen/Dense>

#include "pre.h"
#include "read_config.h"
#include "read_sgrid.h"
#include "create_panel.h"
#include "create_field.h"
#include "map_field.h"
#include "create_wake.h"

using namespace std;
using namespace Eigen;

int pre(bool &symY, double &sRef, double &machInf, double &AoA, Vector3d &vInf,
        Network &bPan, Network &wPan, Field &fPan) {

    //// Initialization
    // Temporary arrays to store grid information
    array<array<double, 3>, 8> box; // Corner points defining domain
    MatrixX3d sGrid ; // Surface grid corner points

    // Path definition
    string configParamName;
    string configParamPath;
    string surfGridName;
    string surfGridPath;

    //// Begin preprocessing
    cout << "*****************************" << endl;
    cout << "*Beginning pre-processing...*" << endl;
    cout << "*****************************" << endl;
    cout << "Files should be located in: /IO/" << endl;
    cout << endl;

    // Paths
    cout << "Enter input parameters file name:" << endl;
    cin >> configParamName;
    configParamPath = "../IO/" + configParamName + ".cfg";
    configParamPath = "/Users/workAC/Documents/Adrien/PhD/Codes/C++/FPMv1/IO/N12.cfg";
    //configParamPath = "C:/Adrien/Work/PhD/Thesis/Codes/C++/FPM/IO/M6.cfg";
    cout << "Enter input surface grid file name:" << endl;
    cin >> surfGridName;
    surfGridPath = "../IO/" +surfGridName + ".pts";
    surfGridPath = "/Users/workAC/Documents/Adrien/PhD/Codes/C++/FPMv1/IO/N12.pts";
    //surfGridPath = "C:/Adrien/Work/PhD/Thesis/Codes/C++/FPM/IO/M6.pts";

    // Read config and grid from files
    read_config(configParamPath, symY, sRef, machInf, AoA, box, fPan);
    read_sgrid(surfGridPath, sGrid, bPan);

    // Set freestream velocity and speed of sound
    vInf(0) = cos(AoA);
    vInf(1) = 0;
    vInf(2) = sin(AoA);

    // Create surface grid
    create_panel(sGrid, bPan);
    create_wake(sGrid, bPan, wPan);

    // Create volume grid and cells mapping
    create_field(box, fPan);
    map_field(sGrid, bPan, fPan);

    //// End preprocessing
    cout << "****************************" << endl;
    cout << "*Pre-processing successful!*" << endl;
    cout << "****************************" << endl;
    return 0;
}