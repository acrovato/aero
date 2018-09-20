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

//// Ray casting and Point In Polygon test
// Check that a given vector does not cross a given (body or wake) panel
//
// References:
// http://geomalgorithms.com/a05-_intersect-1.html
// http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
// https://stackoverflow.com/questions/5188561/signed-angle-between-two-3d-vectors-with-same-origin-within-the-same-plane
//
// Inputs:
// - u0, u1, u2: vector joining center of reference cell to center of adjacent cell
// - w0, w1, w2: vector joining panel center to cell center
// - f0, f1, f2: center of reference cell coordinates
// - n0, n1, n2: panel unit normal
// - v00, v01, v02: panel vertex 1 coordinates
// - v10, v11, v12: panel vertex 2 coordinates
// - v20, v21, v22: panel vertex 3 coordinates
// - v30, v31, v32: panel vertex 4 coordinates
//
// Output:
// - 0 (vector crosses the panel) OR 1 (vector does not cross the panel)

#include <iostream>
#include <Eigen/Dense>
#include "cast_ray_pip.h"

#define TOL 1e-6

using namespace std;
using namespace Eigen;

int cast_ray_pip(double u0, double u1, double u2, double w0, double w1, double w2,
                 double f0, double f1, double f2, double n0, double n1, double n2,
                 double v00, double v01, double v02, double v10, double v11, double v12,
                 double v20, double v21, double v22, double v30, double v31, double v32) {

    double D, N, sI;
    Vector3d a, b, u, w, f, n, i, v0, v1, v2, v3;
    u << u0, u1, u2; // vector joining center of reference cell to center of adjacent cell
    w << w0, w1, w2; // vector joining panel center to cell center
    f << f0, f1, f2; // center of reference cell coordinates
    n << n0, n1, n2; // panel unit normal
    v0 << v00, v01, v02; // panel vertex 1 coordinates
    v1 << v10, v11, v12; // panel vertex 2 coordinates
    v2 << v20, v21, v22; // panel vertex 3 coordinates
    v3 << v30, v31, v32; // panel vertex 4 coordinates

    //// Ray casting

    D = n.dot(u);
    N = -n.dot(w);

    // A. Check if the vector is not parallel to the panel plane
    if (abs(D) < TOL) {
        return 1; // Vector is parallel to, or contained in, the panel plane
    }

    // B. Check if the vector is crossing the panel plane
    sI = N/D; // parameter for intersection
    if (sI < 0 || sI > 1) {
        return 1; // Vector does not reach the panel plane
    }
    // C. Compute the intersection of the vector and the plane
    else {
        i = f + sI * u; // Compute intersection
    }

    //// Point in polygon

    // D. Check that, for each pair of consecutive adjacent vertices, the angle between the lines joining the
    // intersection the vertices is between 0 and PI, measured counterclockwise

    // Vectors from intersection to vertices 1-4
    a = v0 - i;
    b = v3 - i;
    // If (a X b) vector has not the same orientation as panel normal, then angle is between PI and 2PI
    if ((a.cross(b)).dot(n) < -TOL)
        return 1;
    // Vectors from intersection to vertices 4-3
    a = v3 - i;
    b = v2 - i;
    if ((a.cross(b)).dot(n) < -TOL)
        return 1;
    // Vectors from intersection to vertices 3-2
    a = v2 - i;
    b = v1 - i;
    if ((a.cross(b)).dot(n) < -TOL)
        return 1;
    // Vectors from intersection to vertices 2-1
    a = v1 - i;
    b = v0 - i;
    if ((a.cross(b)).dot(n) < -TOL)
        return 1;

    // If program ran till there, then derivative intersects the panel and must be prevented!
    return 0;
}