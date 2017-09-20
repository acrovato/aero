//// Field derivatives mapping
// Identify cells for which spatial derivatives should be prevented so that surface panel crossing does not occur
// Results are stored in boolean matrices: row(f) = [dx-, dx+, dy-, dy+, dz-, dz+]
//
// I/O:
// - sGrid: temporary dynamic array containing body panel vertices
// - numC: list of numerical parameters (structure)
// - bPan: body panels (structure)
// - bPan: wake panels (structure)
// - fPan: field panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "map_derivatives.h"
#include "cast_ray_pip.h"

using namespace std;
using namespace Eigen;

void map_derivatives(MatrixX3d &sGrid, Numerical_CST &numC, Network &bPan, Network &wPan, Field &fPan) {

    //// Begin
    // Variable definition
    int allow = 0;
    double minX, minY, maxY, minZ, maxZ; // minimal bounding box holding the geometry
    Vector3d u, w;
    // Display check
    cout << "Mapping field derivatives... " << flush;
    // Resize field and network matrices
    fPan.fbdMap.resize(fPan.nF,6);
    fPan.fwdMap.resize(fPan.nF,6);
    // Define box enclosing the body
    minX = sGrid.col(0).minCoeff() - numC.TOLB;
    minY = sGrid.col(1).minCoeff() - numC.TOLB;
    maxY = sGrid.col(1).maxCoeff() + numC.TOLB;
    minZ = sGrid.col(2).minCoeff() - numC.TOLB;
    maxZ = sGrid.col(2).maxCoeff() + numC.TOLB;

    //// Map field derivatives
    int f = 0;
    for (int j = 0; j < fPan.nY; j++) {
        for (int k = 0; k < fPan.nZ; k++) {
            for (int i = 0; i < fPan.nX; i++) {

                // If cell is inside the body, derivatives are prevented
                if (!fPan.fMap(f)) {
                    fPan.fbdMap(f, 0) = 0;
                    fPan.fbdMap(f, 1) = 0;
                    fPan.fbdMap(f, 2) = 0;
                    fPan.fbdMap(f, 3) = 0;
                    fPan.fbdMap(f, 4) = 0;
                    fPan.fbdMap(f, 5) = 0;
                    fPan.fwdMap(f, 0) = 0;
                    fPan.fwdMap(f, 1) = 0;
                    fPan.fwdMap(f, 2) = 0;
                    fPan.fwdMap(f, 3) = 0;
                    fPan.fwdMap(f, 4) = 0;
                    fPan.fwdMap(f, 5) = 0;
                }
                else {
                    /// 1. X-backward
                    // Check if cell is not on the domain border
                    if (i == 0) {
                        fPan.fbdMap(f, 0) = 0;
                        fPan.fwdMap(f, 0) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if ((fPan.CG(f,0) < minX)
                             || (fPan.CG(f,1) < minY || fPan.CG(f,1) > maxY)
                             || (fPan.CG(f,2) < minZ || fPan.CG(f,2) > maxZ)) {
                        fPan.fbdMap(f, 0) = 1;
                        fPan.fwdMap(f, 0) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f-1).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 0) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 0) = allow; // Store the result in a boolean matrix
                    }

                    /// 2. X-forward
                    // Check if cell is not on the domain border
                    if (i == fPan.nX - 1) {
                        fPan.fbdMap(f, 1) = 0;
                        fPan.fwdMap(f, 1) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if (((fPan.CG(f,0)+fPan.deltaX) < minX)
                             || (fPan.CG(f,1) < minY || fPan.CG(f,1) > maxY)
                             || (fPan.CG(f,2) < minZ || fPan.CG(f,2) > maxZ)) {
                        fPan.fbdMap(f, 1) = 1;
                        fPan.fwdMap(f, 1) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f+1).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 1) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 1) = allow; // Store the result in a boolean matrix
                    }

                    /// 3. Y-backward
                    // Check if cell is not on the domain border
                    if (j == 0) {
                        fPan.fbdMap(f, 2) = 0;
                        fPan.fwdMap(f, 2) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if ((fPan.CG(f,0) < minX)
                             || (fPan.CG(f,1) < minY || (fPan.CG(f,1)-fPan.deltaY) > maxY)
                             || (fPan.CG(f,2) < minZ || fPan.CG(f,2) > maxZ)) {
                        fPan.fbdMap(f, 2) = 1;
                        fPan.fwdMap(f, 2) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f-fPan.nX*fPan.nZ).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 2) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 2) = allow; // Store the result in a boolean matrix
                    }

                    /// 4. Y-forward
                    // Check if cell is not on the domain border
                    if (j == fPan.nY - 1) {
                        fPan.fbdMap(f, 3) = 0;
                        fPan.fwdMap(f, 3) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if ((fPan.CG(f,0) < minX)
                             || ((fPan.CG(f,1)+fPan.deltaY) < minY || fPan.CG(f,1) > maxY)
                             || (fPan.CG(f,2) < minZ || fPan.CG(f,2) > maxZ)) {
                        fPan.fbdMap(f, 3) = 1;
                        fPan.fwdMap(f, 3) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f+fPan.nX*fPan.nZ).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 3) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 3) = allow; // Store the result in a boolean matrix
                    }

                    /// 5. Z-backward
                    // Check if cell is not on the domain border
                    if (k == 0) {
                        fPan.fbdMap(f, 4) = 0;
                        fPan.fwdMap(f, 4) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if ((fPan.CG(f,0) < minX)
                             || (fPan.CG(f,1) < minY || fPan.CG(f,1) > maxY)
                             || (fPan.CG(f,2) < minZ || (fPan.CG(f,2)-fPan.deltaZ) > maxZ)) {
                        fPan.fbdMap(f, 4) = 1;
                        fPan.fwdMap(f, 4) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f-fPan.nX).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 4) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 4) = allow; // Store the result in a boolean matrix
                    }

                    /// 6. Z-forward
                    // Check if cell is not on the domain border
                    if (k == fPan.nZ - 1) {
                        fPan.fbdMap(f, 5) = 0;
                        fPan.fwdMap(f, 5) = 0;
                    }
                    // Create a box enclosing the geometry, based on cell dimensions to restrict the algorithm complexity
                    else if ((fPan.CG(f,0) < minX)
                             || (fPan.CG(f,1) < minY || fPan.CG(f,1) > maxY)
                             || ((fPan.CG(f,2)+fPan.deltaZ) < minZ || fPan.CG(f,2) > maxZ)) {
                        fPan.fbdMap(f, 5) = 1;
                        fPan.fwdMap(f, 5) = 1;
                    }
                    // Create a vector joining the center of two adjacent cells and check if that vector crosses a panel
                    else {
                        u = fPan.CG.row(f+fPan.nX).transpose() - fPan.CG.row(f).transpose();
                        // Body panel
                        for (int p = 0; p < bPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - bPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 bPan.n(p,0),bPan.n(p,1),bPan.n(p,2),bPan.v0(p,0),bPan.v0(p,1),bPan.v0(p,2),
                                                 bPan.v1(p,0),bPan.v1(p,1),bPan.v1(p,2),bPan.v2(p,0),bPan.v2(p,1),bPan.v2(p,2),
                                                 bPan.v3(p,0),bPan.v3(p,1),bPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fbdMap(f, 5) = allow; // Store the result in a boolean matrix
                        // Wake panel
                        for (int p = 0; p < wPan.nP; p++) {
                            w = fPan.CG.row(f).transpose() - wPan.CG.row(p).transpose();

                            allow = cast_ray_pip(u(0),u(1),u(2),w(0),w(1),w(2),fPan.CG(f,0),fPan.CG(f,1),fPan.CG(f,2),
                                                 wPan.n(p,0),wPan.n(p,1),wPan.n(p,2),wPan.v0(p,0),wPan.v0(p,1),wPan.v0(p,2),
                                                 wPan.v1(p,0),wPan.v1(p,1),wPan.v1(p,2),wPan.v2(p,0),wPan.v2(p,1),wPan.v2(p,2),
                                                 wPan.v3(p,0),wPan.v3(p,1),wPan.v3(p,2));
                            if (!allow)
                                break;
                        }
                        fPan.fwdMap(f, 5) = allow; // Store the result in a boolean matrix
                    }
                }
                f++;
            }
        }
    }

    //// Control display
    cout << "Done!" << endl;
    cout << "Allowed derivatives (body) map: " << fPan.fbdMap.rows() << 'X' << fPan.fbdMap.cols() << endl;
    //for (int i = 0; i < fPan.nF; i++)
    //    cout << fPan.fbdMap.row(i) << endl;
    //cout << endl;
    cout << "Allowed derivatives (wake) map: " << fPan.fwdMap.rows() << 'X' << fPan.fwdMap.cols() << endl;
    //for (int i = 0; i < fPan.nF; i++)
    //    cout << fPan.fwdMap.row(i) << endl;
    cout << endl;
}