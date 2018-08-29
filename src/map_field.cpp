//// Field mapping
// Identify each cell as an interior or exterior (FLAG 0, 1). Also create vector of indexes associated to
// field mapping type.
// Crossing number method is used as point in polygon algorithm.
//
// References: http://geomalgorithms.com/a03-_inclusion.html
//
// I/O:
// - sGrid: temporary dynamic array containing body panel vertices
// - numC: list of numerical parameters (structure)
// - bPan: body panels (structure)
// - fPan: field panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "map_field.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void map_field(MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Field &fPan) {

    //// Begin
    // Temporary variables
    double minX, maxX, minY, maxY, minZ, maxZ; // minimal bounding box holding the geometry
    int iE = 0, iI = 0; // indexes
    MatrixX2d panVert;
    panVert.resize(bPan.nC,2);

    // Display check
    cout << "Mapping field cells... " << flush;

    // Resize matrices
    fPan.fMap.resize(fPan.nF);

    // Compute minimum size box including full geometry
    minX = sGrid.col(0).minCoeff() - numC.TOLB;
    maxX = sGrid.col(0).maxCoeff() + numC.TOLB;
    minY = sGrid.col(1).minCoeff();
    maxY = sGrid.col(1).maxCoeff();
    minZ = sGrid.col(2).minCoeff() - numC.TOLB;
    maxZ = sGrid.col(2).maxCoeff() + numC.TOLB;

    //// Map field cells
    for (int f = 0; f < fPan.nF; f++) {
        lbl_loop0:
        // Bounding box
        if (fPan.CG(f,0) < minX || fPan.CG(f,0) > maxX
            || fPan.CG(f,1) < minY || fPan.CG(f,1) > maxY
            || fPan.CG(f,2) < minZ || fPan.CG(f,2) > maxZ) {
            fPan.fMap(f) = 1;
            fPan.nE++;
        }
        // 2D PIP algorithm
        else {
            // Find spanwise station corresponding to y-coordinate of field cell
            for (int s = 0; s < bPan.nS_; s++) {
                if (fPan.CG(f,1) <= sGrid((s+1)*bPan.nC,1)) {
                    // Interpolate
                    double a = (sGrid((s+1)*bPan.nC,1) - fPan.CG(f,1)) / (sGrid((s+1)*bPan.nC,1) - sGrid(s*bPan.nC,1));
                    double b = (fPan.CG(f,1) - sGrid(s*bPan.nC,1)) / (sGrid((s+1)*bPan.nC,1) - sGrid(s*bPan.nC,1));
                    // Store interpolated airfoil panel vertices
                    for (int i = 0; i < bPan.nC; i++) {
                        panVert(i,0) = a*sGrid(s*bPan.nC+i,0) + b*sGrid((s+1)*bPan.nC+i,0);
                        panVert(i,1) = a*sGrid(s*bPan.nC+i,2) + b*sGrid((s+1)*bPan.nC+i,2);
                    }
                    // Crossing number method
                    int nInter = 0;
                    for (int i = 0; i < bPan.nC_; i++) {
                        // Check if point is not on a vertex
                        if (fPan.CG(f,0) == panVert(i,0) && fPan.CG(f,2) == panVert(i,1)) {
                            fPan.fMap(f) = 0;
                            fPan.nI++;
                            f++;
                            goto lbl_loop0; // *might wanna use goto here!!
                        }
                        // Exclude top endpoint of segment to avoid double crossings and horizontal edges
                        if ((fPan.CG(f,2) >= panVert(i,1) && fPan.CG(f,2) < panVert(i+1,1))
                            || (fPan.CG(f,2) < panVert(i,1) && fPan.CG(f,2) >= panVert(i+1,1))) {
                            double sI = (fPan.CG(f,2) - panVert(i,1)) / (panVert(i+1,1) - panVert(i,1)); // Compute y-intersection
                            if (fPan.CG(f,0) <  panVert(i,0) + sI * (panVert(i+1,0) - panVert(i,0))) { // If point is left of x-intersect
                                nInter++; // Then intersection is valid
                            }
                        }
                    }
                    // Even number of intersections, external point
                    if (nInter%2 == 0) {
                        fPan.fMap(f) = 1;
                        fPan.nE++;
                        break;
                    }
                    // Odd number of intersections, internal point
                    else {
                        fPan.fMap(f) = 0;
                        fPan.nI++;
                        break;
                    }
                }
                else
                    continue;
            }
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

    //// Control display
    cout << "Done!" << endl;
    cout << "Field map: " << fPan.fMap.size() << endl;
    #ifdef VERBOSE
        for (int j = 0; j < fPan.nY; ++j) {
            cout << endl;
            for (int k = 0; k < fPan.nZ; ++k) {
                cout << endl;
                for (int i = 0; i < fPan.nX; ++i)
                    cout << fPan.fMap(i + k * fPan.nX + j * fPan.nX * fPan.nZ) << ' ';
            }
        }
        cout << endl;
    #endif
    cout << "Exterior cells: " << fPan.eIdx.rows() << endl;
    cout << "Interior cells: " << fPan.iIdx.rows() << endl;
    cout << endl;
}