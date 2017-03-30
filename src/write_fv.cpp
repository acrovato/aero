//// Field quantities writing
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
    cout << "Writing field variables file in 'fq.dat'... ";

    // Set path
    outPath += "fq.dat";

    //// Write to file
    ofstream fq;
    fq.open (outPath);

    // General information (header)
    fq << "Field variables file" << endl << endl;
    fq.width(25); fq << left << "Angle of attack: "; fq.width(20); fq << left << alpha*180/3.14159 << endl;
    fq.width(25); fq << left << "Mach: "; fq.width(20); fq << left << Minf << endl;
    fq.width(25); fq << left << "Field cells: "; fq.width(20); fq << left << fPan.nF << endl << endl;

    // Labeling
    fq.width(15); fq << right << "x";
    fq.width(15); fq << right << "y";
    fq.width(15); fq << right << "z";
    fq.width(15); fq << right << "U_x";
    fq.width(15); fq << right << "U_y";
    fq.width(15); fq << right << "U_z";
    fq.width(15); fq << right << "Mach";
    fq.width(15); fq << right << "Sigma" << endl;

    // Write field quantities
    for (int i = 0; i < fPan.nF; ++i) {
        fq.width(15); fq << right << fPan.CG(i,0);
        fq.width(15); fq << right << fPan.CG(i,1);
        fq.width(15); fq << right << fPan.CG(i,2);
        fq.width(15); fq << right << fPan.U(i,0);
        fq.width(15); fq << right << fPan.U(i,1);
        fq.width(15); fq << right << fPan.U(i,2);
        fq.width(15); fq << right << fPan.M(i);
        fq.width(15); fq << right << fPan.sigma(i) << endl;
    }

    // Close file
    fq.close();
    cout << "Done!" << endl << endl;
}
