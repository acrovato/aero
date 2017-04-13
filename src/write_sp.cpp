//// Surface pressure writing
// Write panel pressure surface in external .dat file
//
// I/O:
// - outPath: path to write file
// - sRef: reference surface of the full wing
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - cL: lift coefficient
// - bPan: (network of) body panels
// - cD: drag coefficient

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "write_sp.h"

#define  PI 3.14159

using namespace std;
using namespace Eigen;

void write_sp(string outPath, double sRef, double alpha, double Minf, Network &bPan, double cL, double cD) {

    //// Begin
    cout << "Writing surface pressure file in 'sp.dat'... ";

    // Temporary variables
    int idx = 0;
    outPath += "sp.dat";

    //// Write to file
    ofstream sp;
    sp.open (outPath);

    // General information (header)
    sp << "Surface pressure file" << endl << endl;
    sp.width(25); sp << left << "Angle of attack: "; sp.width(20); sp << left << alpha*180/PI << endl;
    sp.width(25); sp << left << "Mach: "; sp.width(20); sp << left << Minf << endl;
    sp.width(25); sp << left << "Panels: "; sp.width(20); sp << left << bPan.nP << endl;
    sp.width(25); sp << left << "Full reference surface: "; sp.width(20); sp << left << sRef << endl;
    sp.width(25); sp << left << "Lift coefficient: "; sp.width(20); sp << left << cL << endl;
    sp.width(25); sp << left << "Drag coefficient: "; sp.width(20); sp << left << cD << endl;
    sp.width(25); sp << left << "Lift to drag ratio: "; sp.width(20); sp << left << cL/cD << endl << endl;

    // Panel pressure coefficient
    for (int j = 0; j < bPan.nS_; ++j) {
        // Labeling
        sp << "*** Stat" << j << " - y = "<< bPan.CG(idx + 1,1) << " ***" << endl;
        sp.width(15); sp << right << "x";
        sp.width(15); sp << right << "z";
        sp.width(15); sp << right << "cp" << endl;
        for (int i = 0; i < bPan.nC_; ++i) {
            idx = i + j * bPan.nC_;
            sp.width(15); sp << right << bPan.CG(idx,0);
            sp.width(15); sp << right << bPan.CG(idx,2);
            sp.width(15); sp << right << bPan.cP(idx) << endl;
        }
        sp << endl;
    }

    // Close file
    sp.close();
    cout << "Done!" << endl;
}