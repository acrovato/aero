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

//// Interpolation from center to vertices
// Interpolate linearly surface singularities from adjacent panel centers to panel vertices. If panel is on an edge of a
// network, linear extrapolation is used instead
//
// Inputs:
// - idG: current panel global index
// - idC: current panel chordwise index
// - idS: current panel spanwise index
// - bPan: body panels (structure)
//
// Output:
// - vis: Matrix of interpolate singularities (row = vertex number; col 0 = doublet, col 1 = source)

#include <iostream>
#include <Eigen/Dense>
#include "interp_ctv.h"
#include <interp.h>

#define NDIM 3
#define NV 4

using namespace std;
using namespace Eigen;

MatrixXd interp_ctv(int idG, int idC, int idS, Network &bPan) {

    //// Initialization
    // Temporary
    MatrixXd v; // vertices coordinates
    v.resize(NV, NDIM);
    MatrixXi c; // centers (corner of interpolation)
    c.resize(NV, 4);
    MatrixXd vis; // interpolated singularities on vertices
    vis.resize(NV, 2);

    // Store vertices coordinates
    v.row(0) = bPan.v0.row(idG);
    v.row(1) = bPan.v1.row(idG);
    v.row(2) = bPan.v2.row(idG);
    v.row(3) = bPan.v3.row(idG);

    //// Compute panel center indices
    // If panel is on the last row/column, linear extrapolation is used
    if (idS == 0) {
        if (idC == 0) {
            c(0,0) = idG;
            c(0,1) = idG + 1;
            c(0,2) = idG + 1 + bPan.nC_;
            c(0,3) = idG + bPan.nC_;
            c(1,0) = idG;
            c(1,1) = idG + 1;
            c(1,2) = idG + 1 + bPan.nC_;
            c(1,3) = idG + bPan.nC_;
            c(2,0) = idG;
            c(2,1) = idG + 1;
            c(2,2) = idG + 1 + bPan.nC_;
            c(2,3) = idG + bPan.nC_;
            c(3,0) = idG;
            c(3,1) = idG + 1;
            c(3,2) = idG + 1 + bPan.nC_;
            c(3,3) = idG + bPan.nC_;
        }
        else if (idC == bPan.nC_ - 1) {
            c(0,0) = idG - 1;
            c(0,1) = idG;
            c(0,2) = idG + bPan.nC_;
            c(0,3) = idG - 1 + bPan.nC_;
            c(1,0) = idG - 1;
            c(1,1) = idG;
            c(1,2) = idG + bPan.nC_;
            c(1,3) = idG - 1 + bPan.nC_;
            c(2,0) = idG - 1;
            c(2,1) = idG;
            c(2,2) = idG + bPan.nC_;
            c(2,3) = idG - 1 + bPan.nC_;
            c(3,0) = idG - 1;
            c(3,1) = idG;
            c(3,2) = idG + bPan.nC_;
            c(3,3) = idG - 1 + bPan.nC_;
        }
        else {
            c(0,0) = idG - 1;
            c(0,1) = idG;
            c(0,2) = idG + bPan.nC_;
            c(0,3) = idG - 1 + bPan.nC_;
            c(1,0) = idG;
            c(1,1) = idG + 1;
            c(1,2) = idG + 1 + bPan.nC_;
            c(1,3) = idG + bPan.nC_;
            c(2,0) = idG;
            c(2,1) = idG + 1;
            c(2,2) = idG + 1 + bPan.nC_;
            c(2,3) = idG + bPan.nC_;
            c(3,0) = idG - 1;
            c(3,1) = idG;
            c(3,2) = idG + bPan.nC_;
            c(3,3) = idG - 1 + bPan.nC_;
        }
    }
    else if (idS == bPan.nS_ - 1) {
        if (idC == 0) {
            c(0,0) = idG - bPan.nC_;
            c(0,1) = idG - bPan.nC_ + 1;
            c(0,2) = idG + 1;
            c(0,3) = idG;
            c(1,0) = idG - bPan.nC_;
            c(1,1) = idG - bPan.nC_ + 1;
            c(1,2) = idG + 1;
            c(1,3) = idG;
            c(2,0) = idG - bPan.nC_;
            c(2,1) = idG - bPan.nC_ + 1;
            c(2,2) = idG + 1;
            c(2,3) = idG;
            c(3,0) = idG - bPan.nC_;
            c(3,1) = idG - bPan.nC_ + 1;
            c(3,2) = idG + 1;
            c(3,3) = idG;
        }
        else if (idC == bPan.nC_ - 1) {
            c(0,0) = idG - 1 - bPan.nC_;
            c(0,1) = idG - bPan.nC_;
            c(0,2) = idG;
            c(0,3) = idG - 1;
            c(1,0) = idG - 1 - bPan.nC_;
            c(1,1) = idG - bPan.nC_;
            c(1,2) = idG;
            c(1,3) = idG - 1;
            c(2,0) = idG - 1 - bPan.nC_;
            c(2,1) = idG - bPan.nC_;
            c(2,2) = idG;
            c(2,3) = idG - 1;
            c(3,0) = idG - 1 - bPan.nC_;
            c(3,1) = idG - bPan.nC_;
            c(3,2) = idG;
            c(3,3) = idG - 1;
        }
        else {
            c(0,0) = idG - 1 - bPan.nC_;
            c(0,1) = idG - bPan.nC_;
            c(0,2) = idG;
            c(0,3) = idG - 1;
            c(1,0) = idG - bPan.nC_;
            c(1,1) = idG - bPan.nC_ + 1;
            c(1,2) = idG + 1;
            c(1,3) = idG;
            c(2,0) = idG - bPan.nC_;
            c(2,1) = idG - bPan.nC_ + 1;
            c(2,2) = idG + 1;
            c(2,3) = idG;
            c(3,0) = idG - 1 - bPan.nC_;
            c(3,1) = idG - bPan.nC_;
            c(3,2) = idG;
            c(3,3) = idG - 1;
        }
    }
    else {
        c(0,0) = idG - 1 - bPan.nC_;
        c(0,1) = idG - bPan.nC_;
        c(0,2) = idG;
        c(0,3) = idG - 1;
        c(1,0) = idG - bPan.nC_;
        c(1,1) = idG - bPan.nC_ + 1;
        c(1,2) = idG + 1;
        c(1,3) = idG;
        c(2,0) = idG;
        c(2,1) = idG + 1;
        c(2,2) = idG + 1 + bPan.nC_;
        c(2,3) = idG + bPan.nC_;
        c(3,0) = idG - 1;
        c(3,1) = idG;
        c(3,2) = idG + bPan.nC_;
        c(3,3) = idG - 1 + bPan.nC_;
    }

    //// Interpolate each vertex
    for (int l = 0; l < NV; l++) {
        vis(l,0) = interp(bPan.CG(c(l,0),0), bPan.CG(c(l,0),1), bPan.CG(c(l,0),2),
                          bPan.CG(c(l,1),0), bPan.CG(c(l,1),1), bPan.CG(c(l,1),2),
                          bPan.CG(c(l,2),0), bPan.CG(c(l,2),1), bPan.CG(c(l,2),2),
                          bPan.CG(c(l,3),0), bPan.CG(c(l,3),1), bPan.CG(c(l,3),2),
                          bPan.mu(c(l,0)), bPan.mu(c(l,1)), bPan.mu(c(l,2)), bPan.mu(c(l,3)),
                          v(l,0), v(l,1), v(l,2));

        vis(l,1) = interp(bPan.CG(c(l,0),0), bPan.CG(c(l,0),1), bPan.CG(c(l,0),2),
                          bPan.CG(c(l,1),0), bPan.CG(c(l,1),1), bPan.CG(c(l,1),2),
                          bPan.CG(c(l,2),0), bPan.CG(c(l,2),1), bPan.CG(c(l,2),2),
                          bPan.CG(c(l,3),0), bPan.CG(c(l,3),1), bPan.CG(c(l,3),2),
                          bPan.tau(c(l,0)), bPan.tau(c(l,1)), bPan.tau(c(l,2)), bPan.tau(c(l,3)),
                          v(l,0), v(l,1), v(l,2));
    }
    return vis;
}