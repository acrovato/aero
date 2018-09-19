//// Field variables writing
// Write field variable in external .dat and .pos files
//
// I/O:
// - outPath: path to write file
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - fPan: (network of) field panels

/* Copyright (C) 2018 Adrien Crovato */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "write_fv.h"

using namespace std;
using namespace Eigen;

void write_fv(string outPath, double alpha, double Minf, Field &fPan) {

    //// Write to .dat file
    cout << "Writing field variables file in 'fv.dat'... " << flush;
    string outPathDat = outPath + "fv.dat";
    ofstream fv;
    fv.open(outPathDat);

    // General information (header)
    fv << "Field variables file" << endl << endl;
    fv.width(25); fv << left << "Angle of attack: "; fv.width(20); fv << left << alpha*180/3.14159 << endl;
    fv.width(25); fv << left << "Mach: "; fv.width(20); fv << left << Minf << endl;
    fv.width(25); fv << left << "Field cells: "; fv.width(20); fv << left << fPan.nF << endl << endl;

    // Labeling
    fv.width(15); fv << right << "x";
    fv.width(15); fv << right << "y";
    fv.width(15); fv << right << "z";
    fv.width(15); fv << right << "U_x";
    fv.width(15); fv << right << "U_y";
    fv.width(15); fv << right << "U_z";
    fv.width(15); fv << right << "rho";
    fv.width(15); fv << right << "rho_X";
    fv.width(15); fv << right << "rho_Z";
    fv.width(15); fv << right << "Mach";
    fv.width(15); fv << right << "Sigma";
    fv.width(15); fv << right << "Epsilon" << endl;

    // Write field quantities
    for (int i = 0; i < fPan.nF; ++i) {
        fv.width(15); fv << right << fPan.CG(i,0);
        fv.width(15); fv << right << fPan.CG(i,1);
        fv.width(15); fv << right << fPan.CG(i,2);
        fv.width(15); fv << right << fPan.U(i,0);
        fv.width(15); fv << right << fPan.U(i,1);
        fv.width(15); fv << right << fPan.U(i,2);
        fv.width(15); fv << right << fPan.rho(i);
        fv.width(15); fv << right << fPan.dRho(i,0);
        fv.width(15); fv << right << fPan.dRho(i,2);
        fv.width(15); fv << right << fPan.M(i);
        fv.width(15); fv << right << fPan.sigma(i);
        fv.width(15); fv << right << fPan.epsilon(i) << endl;
    }

    // Close file
    fv.close();
    cout << "Done!" << endl;

    //// Write to gmsh .pos file
    cout << "Writing field variables file in 'M.pos'... " << flush;
    string outPathPos = outPath + "M.pos";
    fv.open(outPathPos);

    // Gmsh header
    fv << "View \"M\" {" << endl;

    // Mach number
    for (int f = 0; f < fPan.nF; ++f) {
        fv << "SH(";
        fv << fPan.vX(f,0) << "," << fPan.vY(f,0) << "," << fPan.vZ(f,0) << ',';
        fv << fPan.vX(f,1) << "," << fPan.vY(f,0) << "," << fPan.vZ(f,0) << ',';
        fv << fPan.vX(f,1) << "," << fPan.vY(f,1) << "," << fPan.vZ(f,0) << ',';
        fv << fPan.vX(f,0) << "," << fPan.vY(f,1) << "," << fPan.vZ(f,0) << ',';
        fv << fPan.vX(f,0) << "," << fPan.vY(f,0) << "," << fPan.vZ(f,1) << ',';
        fv << fPan.vX(f,1) << "," << fPan.vY(f,0) << "," << fPan.vZ(f,1) << ',';
        fv << fPan.vX(f,1) << "," << fPan.vY(f,1) << "," << fPan.vZ(f,1) << ',';
        fv << fPan.vX(f,0) << "," << fPan.vY(f,1) << "," << fPan.vZ(f,1) << "){";
        for (int i = 0; i < 7; ++i)
            fv << fPan.M(f) << ",";
        fv << fPan.M(f) << "};" << endl;
    }

    // Gmsh footer
    fv << "};" << endl << endl;

    // Close file
    fv.close();
    cout << "Done!" << endl << endl;
}
