//// Post-processing
// Write surface and field quantities to file
//
// Inputs:
// - sRef: reference surface of the full wing
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - bPan: (network of) body panels
// - fPan: (network of) field panels
// - cL: lift coefficient
// - cD: drag coefficient
//
// Outputs:
// - 0, if function succeeded
// - 1, if function failed

#include <iostream>
#include <Eigen/Dense>
#include "post.h"
#include "write_sp.h"
#include "write_fv.h"

using namespace std;
using namespace Eigen;

int post(double sRef, double alpha, double Minf, Network &bPan, Field &fPan, double cL, double cD) {

    //// Begin post-processing
    cout << "******************************" << endl;
    cout << "*Beginning post-processing...*" << endl;
    cout << "******************************" << endl << endl;

    // Paths
    string outPath;
    outPath = "../IO/";

    //// Write to file
    write_sp(outPath, sRef, alpha, Minf, bPan, cL, cD);
    write_fv(outPath, alpha, Minf, fPan);

    //// End post-processing
    cout << "*****************************" << endl;
    cout << "*Post-processing successful!*" << endl;
    cout << "*****************************" << endl;
    return 0;
}