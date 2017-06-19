//// Field creation
// Create field cells from data contained into box. Divide the initial box (whose coordinates are provided as input)
// by the number of cells (also provided as input). Compute center and corner points of field cells.
//
// I/O:
// - box: temporary array containg vertices of domain
// - fPan: field panels (structure)

#include <array>
#include <iostream>
#include <Eigen/Dense>
#include "create_field.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void create_field(array<array<double, 3>, 8> &box, Field &fPan) {

    // Temporary variables
    VectorXd X(fPan.nX,1), Y(fPan.nY,1), Z(fPan.nZ,1);
    VectorXd X0(fPan.nX,1), Y0(fPan.nY,1), Z0(fPan.nZ,1);
    VectorXd X1(fPan.nX,1), Y1(fPan.nY,1), Z1(fPan.nZ,1);
    int idx = 0;

    //// Begin
    cout << "Creating field cells... " << flush;

    // Resize field matrices
    fPan.CG.resize(fPan.nF, NDIM);
    fPan.vX.resize(fPan.nF, 2);
    fPan.vY.resize(fPan.nF, 2);
    fPan.vZ.resize(fPan.nF, 2);

    // Divide the domain
    fPan.deltaX = abs(box[1][0] - box[0][0]) / fPan.nX;
    for (int i = 0; i < fPan.nX; ++i) {
        X(i) = box[0][0] + (i + 0.5) * fPan.deltaX;
        X0(i) = box[0][0] + i * fPan.deltaX;
        X1(i) = box[0][0] + (i + 1) * fPan.deltaX;
    }
    fPan.deltaY = (box[4][1] - box[0][1]) / fPan.nY;
    for (int j = 0; j < fPan.nY; ++j) {
        Y(j) = box[0][1] + (j + 0.5) * fPan.deltaY;
        Y0(j) = box[0][1] + j * fPan.deltaY;
        Y1(j) = box[0][1] + (j + 1) * fPan.deltaY;
    }
    fPan.deltaZ = (box[2][2] - box[0][2]) / fPan.nZ;
    for (int k = 0; k < fPan.nZ; ++k) {
        Z(k) = box[0][2] + (k + 0.5) * fPan.deltaZ;
        Z0(k) = box[0][2] + k * fPan.deltaZ;
        Z1(k) = box[0][2] + (k + 1) * fPan.deltaZ;
    }

    // Center points & corner points
    for (int j = 0; j < fPan.nY; ++j) {
        for (int k = 0; k < fPan.nZ; ++k) {
            for (int i = 0; i < fPan.nX; ++i) {
                idx = i + k * fPan.nX + j * fPan.nX * fPan.nZ;
                fPan.CG(idx,0) = X(i);
                fPan.CG(idx,1) = Y(j);
                fPan.CG(idx,2) = Z(k);
                fPan.vX(idx,0) = X0(i);
                fPan.vX(idx,1) = X1(i);
                fPan.vY(idx,0) = Y0(j);
                fPan.vY(idx,1) = Y1(j);
                fPan.vZ(idx,0) = Z0(k);
                fPan.vZ(idx,1) = Z1(k);
            }
        }
    }

    //// Control display
    cout << "Done!" << endl;
    cout << "Field cell center points: " << fPan.CG.rows() << 'X' << fPan.CG.cols() << endl;
    //for (int i = 0; i < div; ++i)
    //    cout << i << ' ' << fldPts(i,0) << ' ' << fldPts(i,1) << ' ' << fldPts(i,2) << endl;
    cout << "Field cell corner points: " << 8 * fPan.vX.rows() << 'X' << fPan.vX.cols() << endl;
    //for (int i = 0; i < div; ++i) {
    //    cout << i << '\t' << fldXCrn(i,0) << ' ' << fldYCrn(i,0) << ' ' << fldZCrn(i,0) << endl;
    //    cout << '\t' << fldXCrn(i,1) << ' ' << fldYCrn(i,0) << ' ' << fldZCrn(i,0) << endl;
    //    cout << '\t' << fldXCrn(i,0) << ' ' << fldYCrn(i,0) << ' ' << fldZCrn(i,1) << endl;
    //    cout << '\t' << fldXCrn(i,1) << ' ' << fldYCrn(i,0) << ' ' << fldZCrn(i,1) << endl;
    //    cout << '\t' << fldXCrn(i,0) << ' ' << fldYCrn(i,1) << ' ' << fldZCrn(i,0) << endl;
    //    cout << '\t' << fldXCrn(i,1) << ' ' << fldYCrn(i,1) << ' ' << fldZCrn(i,0) << endl;
    //    cout << '\t' << fldXCrn(i,0) << ' ' << fldYCrn(i,1) << ' ' << fldZCrn(i,1) << endl;
    //    cout << '\t' << fldXCrn(i,1) << ' ' << fldYCrn(i,1) << ' ' << fldZCrn(i,1) << endl;
    //}
    cout << endl;
}