/* 
// Copyright 2018 University of Liege
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors:
// - Adrien Crovato
*/

//// Wake panels creation
// Compute wake panel collocation point, surface, vertices, normal, longitudinal, transverse and perpendicular
// vectors from data contained into sGrid
// Wake panels extend to a constant x coordinate, based on the geometry
//
// I/O:
// - sGrid: temporary dynamic array containing body panel vertices
// - bPan: body panels (structure)
// - wPan: wake panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "create_wake.h"

#define  NDIM 3

using namespace std;
using  namespace Eigen;

void create_wake(MatrixX3d &sGrid, Network &bPan, Network &wPan) {

    // Temporary variables
    double norm = 0;
    Vector3d v1(NDIM), v2(NDIM);
    double dwnstrm;
    int c2, c3;

    //// Begin
    cout << "Creating wake panels... " << flush;

    // Set number of panels (currently only handles predefined wakes)
    wPan.nC = 2;
    wPan.nC_ = 1;
    wPan.nS = bPan.nS;
    wPan.nS_ = bPan.nS_;
    wPan.nP = wPan.nS_;

    // Resize network matrices
    wPan.CG.resize(wPan.nP, NDIM);
    wPan.v0.resize(wPan.nP, NDIM);
    wPan.v1.resize(wPan.nP, NDIM);
    wPan.v2.resize(wPan.nP, NDIM);
    wPan.v3.resize(wPan.nP, NDIM);
    wPan.S.resize(wPan.nP, 1);
    wPan.n.resize(wPan.nP, NDIM);
    wPan.l.resize(wPan.nP, NDIM);
    wPan.t.resize(wPan.nP, NDIM);
    wPan.p.resize(wPan.nP, NDIM);

    // Downstream coordinate of wake panels = tip x-coord+ 5*chord
    dwnstrm = sGrid(sGrid.rows()-1,0) + 5*sGrid(0,0);

    // Collocations points, corners and vectors
    for (int j = 0; j < wPan.nP; ++j) {

        c2 = j * bPan.nC;
        c3 = (j+1) * bPan.nC;

        wPan.v0(j,0) = dwnstrm; // pts[i][k]
        wPan.v1(j,0) = sGrid(c2,0); // pts[i+1][k]
        wPan.v2(j,0) = sGrid(c3,0); // pts[i+1][k+1]
        wPan.v3(j,0) = dwnstrm; // pts[i][k+1]
        wPan.v0(j,1) = sGrid(c2,1);
        wPan.v1(j,1) = sGrid(c2,1);
        wPan.v2(j,1) = sGrid(c3,1);
        wPan.v3(j,1) = sGrid(c3,1);
        wPan.v0(j,2) = sGrid(c2,2);
        wPan.v1(j,2) = sGrid(c2,2);
        wPan.v2(j,2) = sGrid(c3,2);
        wPan.v3(j,2) = sGrid(c3,2);

        wPan.CG(j,0) = (wPan.v0(j,0) + wPan.v1(j,0) + wPan.v2(j,0) + wPan.v3(j,0)) / 4;
        wPan.CG(j,1) = (wPan.v1(j,1) + wPan.v2(j,1)) / 2;
        wPan.CG(j,2) = (wPan.v1(j,2) + wPan.v2(j,2)) / 2;

        wPan.l(j,0) = (wPan.v0(j,0) + wPan.v3(j,0) - wPan.v1(j,0) - wPan.v2(j,0)) / 4;
        wPan.l(j,1) = 0;
        wPan.l(j,2) = 0;
        wPan.l.row(j) /= wPan.l.row(j).norm();

        wPan.t(j,0) = (wPan.v2(j,0) - wPan.v1(j,0)) / 2;
        wPan.t(j,1) = (wPan.v2(j,1) - wPan.v1(j,1)) / 2;
        wPan.t(j,2) = (wPan.v2(j,2) - wPan.v1(j,2)) / 2;
        wPan.t.row(j) /= wPan.t.row(j).norm();
    }
    // Surfaces and normals
    for (int j = 0; j < wPan.nP; ++j) {
        v1(0) = wPan.v2(j,0) - wPan.v0(j,0);
        v2(0) = wPan.v1(j,0) - wPan.v3(j,0);
        v1(1) = wPan.v2(j,1) - wPan.v0(j,1);
        v2(1) = wPan.v1(j,1) - wPan.v3(j,1);
        v1(2) = wPan.v2(j,2) - wPan.v0(j,2);
        v2(2) = wPan.v1(j,2) - wPan.v3(j,2);

        wPan.n.row(j) = v1.cross(v2);
        norm = wPan.n.row(j).norm();

        wPan.S(j) = norm/2;
        wPan.n.row(j) /= norm;
    }
    // Perpendicular vector
    for (int j = 0; j < wPan.nP; ++j)
        wPan.p.row(j) = wPan.n.row(j).cross(wPan.l.row(j));

    //// Control display
    cout << "Done!" << endl;
    cout << "Wake panels extend to " << dwnstrm << " in the x direction" << endl;
    #ifdef VERBOSE
        cout << "Collocation points: " << wPan.CG.rows() << 'X' << wPan.CG.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.CG(i,0) << ' ' << wPan.CG(i,1) << ' ' << wPan.CG(i,2) << endl;
        cout << "Panel surfaces: " << wPan.S.rows() << 'X' << wPan.S.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.S(i) << endl;
        cout << "Unit normals: " << wPan.n.rows() << 'X' << wPan.n.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.n(i,0) << ' ' << wPan.n(i,1) << ' ' << wPan.n(i,2) << endl;
        cout << "Unit longitudinal vectors: " << wPan.l.rows() << 'X' << wPan.l.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.l(i,0) << ' ' << wPan.l(i,1) << ' ' << wPan.l(i,2) << endl;
        cout << "Unit transverse vectors: " << wPan.t.rows() << 'X' << wPan.t.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.t(i,0) << ' ' << wPan.t(i,1) << ' ' << wPan.t(i,2) << endl;
        cout << "Unit perpendicular vectors: " << wPan.p.rows() << 'X' << wPan.p.cols() << endl;
        for (int i = 0; i < wPan.nP; ++i)
            cout << i << ' ' << wPan.p(i,0) << ' ' << wPan.p(i,1) << ' ' << wPan.p(i,2) << endl;
    #endif
    cout << endl;
}