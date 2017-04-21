//// Map field
// Identify each cell as an interior or exterior (FLAG 0, 1). Also create vector of indexes associated to
// field mapping type. Identify cells just above/below wake (FLAG 1 in wakeMap)
//
// I/O:
// - sGrid: temporary dynamic array containing body panel vertices
// - bPan: body panels (structure)
// - fPan: field panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "map_field.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void map_field(MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Field &fPan) {

    // Temporary variables
    int idx = 0, idx0 = 0, idx1 = 0; // counters
    double a, b, xL, xL0, zL; // parameters for linear interpolation
    Vector3d fCtr, pt, nrm, v; // field center, panel center, panel normal and fCtr-pt vector
    VectorXd dotPrd, dist; // dot product and distance
    double minX, maxX, minY, maxY, minZ, maxZ; // minimal bounding box of holding the geometry
    int iE = 0, iI = 0, iW = 0; // indexes
    dist.resize(bPan.nP);
    dotPrd.resize(bPan.nP);

    //// Begin
    cout << "Mapping field cells... " << flush;

    // Resize field and network matrices
    fPan.fMap.resize(fPan.nF);
    fPan.wMap.resize(fPan.nF);

    // Compute minimum size box including full geometry
    minX = sGrid.col(0).minCoeff() - numC.TOLB;
    maxX = sGrid.col(0).maxCoeff() + numC.TOLB;
    minY = sGrid.col(1).minCoeff() - numC.TOLB;
    maxY = sGrid.col(1).maxCoeff() + numC.TOLB;
    minZ = sGrid.col(2).minCoeff() - numC.TOLB;
    maxZ = sGrid.col(2).maxCoeff() + numC.TOLB;

    //// Map cells
    for (int f = 0; f < fPan.nF; ++f) {
        fCtr = fPan.CG.row(f).transpose();
        if (((fCtr(0) + fPan.deltaMG) >= minX && (fCtr(0) - fPan.deltaMG) <= maxX)
            && ((fCtr(1) + fPan.deltaMG) >= minY && (fCtr(1) - fPan.deltaMG) <= maxY)
            && ((fCtr(2) + fPan.deltaMG) >= minZ && (fCtr(2) - fPan.deltaMG) <= maxZ)) {
            for (int p = 0; p < bPan.nP; ++p) {
                // tmp variables for compatibility
                pt = bPan.CG.row(p).transpose();
                nrm = bPan.n.row(p).transpose();
                // dot product computation
                v = fCtr - pt;
                dotPrd(p) = v.dot(nrm);
                dist(p) = v.norm();
            }
            // if dot product sign changes, point outside surface
            if (dotPrd.maxCoeff() >= 0 && dotPrd.minCoeff() < 0) {
                // Find spanwise location of field cell to handle point to close to surface
                for (int s = 0; s < bPan.nS_; ++s) {
                    idx0 = s * bPan.nC;
                    idx1 = (s + 1) * bPan.nC;
                    if (fPan.CG(f, 1) <= sGrid(idx1, 1)) {
                        a = (sGrid(idx1, 0) - sGrid(idx0, 0)) / (sGrid(idx1, 1) - sGrid(idx0, 1));
                        b = sGrid(idx0, 0) - a * sGrid(idx0, 1);
                        xL = a * fPan.CG(f, 1) + b;
                        a = (sGrid(idx1 + bPan.nC_/2, 0) - sGrid(idx0 + bPan.nC_/2, 0))
                            / (sGrid(idx1 + bPan.nC_/2, 1) - sGrid(idx0 + bPan.nC_/2, 1));
                        b = sGrid(idx0 + bPan.nC_/2, 0) - a * sGrid(idx0 + bPan.nC_/2, 1);
                        xL0 = a * fPan.CG(f, 1) + b;
                        // if cell is NOT fwd or aft local airfoil
                        if (fPan.CG(f,0) + fPan.deltaMG > (xL0-numC.TOLS) && fPan.CG(f,0) - fPan.deltaMG < (xL+numC.TOLS)) {
                            dist.block(s * bPan.nC_, 0, bPan.nC_, 1).minCoeff(&idx);
                            // if normal distance between field point and closest panel is less than TOLS, interior point
                            if (dotPrd(idx + s*bPan.nC_) - fPan.deltaMG < numC.TOLS) {
                                fPan.fMap(f) = 0;
                                fPan.nI++;
                                break;
                            }
                            else {
                                fPan.fMap(f) = 1;
                                fPan.nE++;
                                break;
                            }
                        }
                        else {
                            fPan.fMap(f) = 1;
                            fPan.nE++;
                            break;
                        }
                    }
                    else
                        continue;
                }
            }
            // if dot product sign < 0, point inside surface, interior point
            else if (dotPrd.maxCoeff() < 0 && dotPrd.minCoeff() < 0) {
                fPan.fMap(f) = 0;
                fPan.nI++;
            }
            else {
                cout << endl << "Field cell " << f << " cannot be mapped!" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else {
            fPan.fMap(f) = 1;
            fPan.nE++;
        }
    }

    // Find cell indices according to mapping type
    fPan.eIdx.resize(fPan.nE);
    fPan.iIdx.resize(fPan.nI);
    for (int f = 0; f < fPan.nF; ++f) {
        if (fPan.fMap(f)) {
            fPan.eIdx(iE) = f;
            iE++;
        }
        else {
            fPan.iIdx(iI) = f;
            iI++;
        }
    }

    //// Wake map
    fPan.nW = 0;
    // Find cells just above/below the wake
    for (int f = 0; f < fPan.nF; ++f) {
        if (fPan.CG(f,0) > minX && fPan.CG(f,1) > minY && fPan.CG(f,1) < maxY) {
            for (int p = 0; p < bPan.nS_; ++p) {
                idx0 = p * bPan.nC;
                idx1 = (p + 1) * bPan.nC;
                if (fPan.CG(f,1) < sGrid(idx1,1)) {
                    a = (sGrid(idx1,0) - sGrid(idx0,0)) / (sGrid(idx1,1) - sGrid(idx0,1));
                    b = sGrid(idx0,0) - a*sGrid(idx0,1);
                    xL = a * fPan.CG(f,1) + b;
                    a = (sGrid(idx1,2) - sGrid(idx0,2)) / (sGrid(idx1,1) - sGrid(idx0,1));
                    b = sGrid(idx0,2) - a*sGrid(idx0,1);
                    zL = a * fPan.CG(f,1) + b;
                    if (fPan.CG(f,0) > xL) {
                        if (abs(fPan.CG(f,2) - zL) < 1.01 * fPan.deltaMG) {
                            fPan.wMap(f) = 1;
                            fPan.nW++;
                            break;
                        }
                        else {
                            fPan.wMap(f) = 0;
                            break;
                        }
                    }
                    else {
                        fPan.wMap(f) = 0;
                        break;
                    }
                }
                else
                    continue;
            }
        }
        else {
            fPan.wMap(f) = 0;
            continue;
        }
    }

    // Find cell indices according to wake type
    fPan.wIdx.resize(fPan.nW);
    for (int f = 0; f < fPan.nF; ++f) {
        if (fPan.wMap(f)) {
            fPan.wIdx(iW) = f;
            iW++;
        }
        else
            continue;
    }

    //// Control display
    cout << "Done!" << endl;
    cout << "Field map: " << fPan.fMap.size() << endl;
    //for (int j = 0; j < nY; ++j) {
    //    cout << endl;
    //    for (int k = 0; k < nZ; ++k) {
    //        cout << endl;
    //        for (int i = 0; i < nX; ++i)
    //            cout << fldMap(i + k * nX + j * nX * nZ) << ' ';
    //    }
    //}
    //cout << endl;
    cout << "Exterior cells: " << fPan.eIdx.rows() << endl;
    cout << "Interior cells: " << fPan.iIdx.rows() << endl;
    cout << "Wake map: " << fPan.wMap.size() << endl;
    //for (int j = 0; j < nY; ++j) {
    //    cout << endl;
    //    for (int k = 0; k < nZ; ++k) {
    //        cout << endl;
    //        for (int i = 0; i < nX; ++i)
    //            cout << wakeMap(i + k * nX + j * nX * nZ) << ' ';
    //    }
    //}
    cout << "Wake cells: " << fPan.wIdx.rows() << endl;
    cout << endl;
}