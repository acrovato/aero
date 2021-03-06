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

//// Body panels creation
// Compute panel collocation point, surface, vertices, normal, longitudinal, transverse and perpendicular unit vectors
// from data contained into sGrid
//
// I/O:
// - sGrid: temporary dynamic array containing body panel vertices
// - bPan: body panels (structure)

#include <iostream>
#include <Eigen/Dense>
#include "create_panel.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void create_panel(MatrixX3d &sGrid, Network &bPan) {

    // Temporary variables
    double norm = 0;
    int panIdx;
    int c1, c2, c3, c4;
    Vector3d v1(NDIM), v2(NDIM);

    //// Begin
    cout << "Creating panels... " << flush;

    // Resizing network matrices
    bPan.CG.resize(bPan.nP, NDIM);
    bPan.v0.resize(bPan.nP, NDIM);
    bPan.v1.resize(bPan.nP, NDIM);
    bPan.v2.resize(bPan.nP, NDIM);
    bPan.v3.resize(bPan.nP, NDIM);
    bPan.S.resize(bPan.nP, 1);
    bPan.n.resize(bPan.nP, NDIM);
    bPan.l.resize(bPan.nP, NDIM);
    bPan.t.resize(bPan.nP, NDIM);
    bPan.p.resize(bPan.nP, NDIM);

    // Compute corner points, collocation points and vectors
    for (int i = 0; i < bPan.nC_; ++i) {
        for (int k = 0; k < bPan.nS_; ++k) {

            panIdx = i + k * (bPan.nC_);
            c1 = i + k * bPan.nC;
            c2 = i + 1 + k * bPan.nC;
            c3 = i + 1 + (k+1) *bPan.nC;
            c4 = i + (k+1) * bPan.nC;

            bPan.v0(panIdx,0) = sGrid(c1,0); // pts[i][k]
            bPan.v1(panIdx,0) = sGrid(c2,0); // pts[i+1][k]
            bPan.v2(panIdx,0) = sGrid(c3,0); // pts[i+1][k+1]
            bPan.v3(panIdx,0) = sGrid(c4,0); // pts[i][k+1]
            bPan.v0(panIdx,1) = sGrid(c1,1);
            bPan.v1(panIdx,1) = sGrid(c2,1);
            bPan.v2(panIdx,1) = sGrid(c3,1);
            bPan.v3(panIdx,1) = sGrid(c4,1);
            bPan.v0(panIdx,2) = sGrid(c1,2);
            bPan.v1(panIdx,2) = sGrid(c2,2);
            bPan.v2(panIdx,2) = sGrid(c3,2);
            bPan.v3(panIdx,2) = sGrid(c4,2);

            bPan.CG(panIdx,0) = (bPan.v0(panIdx,0) + bPan.v1(panIdx,0) + bPan.v2(panIdx,0) + bPan.v3(panIdx,0)) / 4;
            bPan.CG(panIdx,1) = (bPan.v0(panIdx,1) + bPan.v1(panIdx,1) + bPan.v2(panIdx,1) + bPan.v3(panIdx,1)) / 4;
            bPan.CG(panIdx,2) = (bPan.v0(panIdx,2) + bPan.v1(panIdx,2) + bPan.v2(panIdx,2) + bPan.v3(panIdx,2)) / 4;

            bPan.l(panIdx,0) = (bPan.v0(panIdx,0) + bPan.v3(panIdx,0) - bPan.v1(panIdx,0) - bPan.v2(panIdx,0)) / 4;
            bPan.l(panIdx,1) = (bPan.v0(panIdx,1) + bPan.v3(panIdx,1) - bPan.v1(panIdx,1) - bPan.v2(panIdx,1)) / 4;
            bPan.l(panIdx,2) = (bPan.v0(panIdx,2) + bPan.v3(panIdx,2) - bPan.v1(panIdx,2) - bPan.v2(panIdx,2)) / 4;
            bPan.l.row(panIdx) /= bPan.l.row(panIdx).norm();

            bPan.t(panIdx,0) = (bPan.v2(panIdx,0) + bPan.v3(panIdx,0) - bPan.v0(panIdx,0) - bPan.v1(panIdx,0)) / 4;
            bPan.t(panIdx,1) = (bPan.v2(panIdx,1) + bPan.v3(panIdx,1) - bPan.v0(panIdx,1) - bPan.v1(panIdx,1)) / 4;
            bPan.t(panIdx,2) = (bPan.v2(panIdx,2) + bPan.v3(panIdx,2) - bPan.v0(panIdx,2) - bPan.v1(panIdx,2)) / 4;
            bPan.t.row(panIdx) /= bPan.t.row(panIdx).norm();
        }
    }
    // Compute panel normals and surfaces
    for (int i = 0; i < bPan.nC_; ++i) {
        for (int k = 0; k < bPan.nS_; ++k) {
            panIdx = i + k * (bPan.nC_);

            v1(0) = bPan.v2(panIdx,0) - bPan.v0(panIdx,0);
            v2(0) = bPan.v1(panIdx,0) - bPan.v3(panIdx,0);
            v1(1) = bPan.v2(panIdx,1) - bPan.v0(panIdx,1);
            v2(1) = bPan.v1(panIdx,1) - bPan.v3(panIdx,1);
            v1(2) = bPan.v2(panIdx,2) - bPan.v0(panIdx,2);
            v2(2) = bPan.v1(panIdx,2) - bPan.v3(panIdx,2);

            bPan.n.row(panIdx) = v1.cross(v2);
            norm = bPan.n.row(panIdx).norm();

            bPan.S(panIdx) = norm/2;
            bPan.n.row(panIdx) /= norm;
        }
    }
    // Compute panel perpendicular vector
    for (int i = 0; i < bPan.nC_; ++i) {
        for (int k = 0; k < bPan.nS_; ++k) {
            panIdx = i + k * (bPan.nC_);
            bPan.p.row(panIdx) = bPan.n.row(panIdx).cross(bPan.l.row(panIdx));
        }
    }

    //// Control display
    cout << "Done!" << endl;
    #ifdef VERBOSE
        cout << "Collocation points: " << bPan.CG.rows() << 'X' << bPan.CG.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.CG(i,0) << ' ' << bPan.CG(i,1) << ' ' << bPan.CG(i,2) << endl;
        cout << "Panel surfaces: " << bPan.S.rows() << 'X' << bPan.S.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.S(i) << endl;
        cout << "Unit normals: " << bPan.n.rows() << 'X' << bPan.n.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.n(i,0) << ' ' << bPan.n(i,1) << ' ' << bPan.n(i,2) << endl;
        cout << "Unit longitudinal vectors: " << bPan.l.rows() << 'X' << bPan.l.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.l(i,0) << ' ' << bPan.l(i,1) << ' ' << bPan.l(i,2) << endl;
        cout << "Unit transverse vectors: " << bPan.t.rows() << 'X' << bPan.t.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.t(i,0) << ' ' << bPan.t(i,1) << ' ' << bPan.t(i,2) << endl;
        cout << "Unit perpendicular vectors: " << bPan.p.rows() << 'X' << bPan.p.cols() << endl;
        for (int i = 0; i < bPan.nP; ++i)
            cout << i << ' ' << bPan.p(i,0) << ' ' << bPan.p(i,1) << ' ' << bPan.p(i,2) << endl;
    #endif
    cout << endl;
}