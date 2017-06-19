//// Field variables writing
// Write field variable in external .dat file
//
// I/O:
// - outPath: path to write file
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - fPan: (network of) field panels

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "write_fv.h"

using namespace std;
using namespace Eigen;

void write_fv(string outPath, double alpha, double Minf, Field &fPan) {

    //// Begin
    cout << "Writing field variables file in 'fv.dat'... ";

    // Set path
    outPath += "fv.dat";

    //// Write to file
    ofstream fv;
    fv.open (outPath);

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
    cout << "Done!" << endl << endl;
}
