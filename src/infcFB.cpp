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

//// Field to body velocity influence coefficient computation
// Compute field source velocity influence coefficient between a field panel and a point
//
// Reference: Chu, L., Yates, E., & Kandil, O. (1989),
// Integral Equation Solution of the Full Potential Equation for Transonic Flows, 27th Aerospace Sciences Meeting.
//
// Inputs:
// - x, y, z: coordinates of target point
// - xC: x coordinates of influencing field cell vertices
// - yC: y coordinates of influencing field cell vertices
// - zC: z coordinates of influencing field cell vertices
//
// Output:
// - coeff: array of AIC ([0] = u, [1] = v, [2] = w)

#include <iostream>
#include <array>
#include <cmath>
#include "infcFB.h"

#define PI 3.14159

using namespace std;

array<double,3> infcFB(double x, double y, double z, array<double,2> &xC, array<double,2> &yC, array<double,2> &zC) {

    // Temporary variables
    array<double,3> coeff;
    double R, A, B, C, t0, t1, t2; // Coefficients
    coeff[0] = 0;
    coeff[1] = 0;
    coeff[2] = 0;

    // AIC
    for (int i = 1; i <= 2; ++i) {
        for (int j = 1; j <= 2; ++j) {
            for (int k = 1; k <= 2; ++k) {
                A = (x-xC[i-1]);
                B = (y-yC[j-1]);
                C = (z-zC[k-1]);
                R = sqrt(A*A + B*B + C*C);
                t0 = B/2 * log((R+C)/(R-C)) + C/2 * log((R+B)/(R-B)) - A * atan((B*C)/(A*R));
                t1 = C/2 * log((R+A)/(R-A)) + A/2 * log((R+C)/(R-C)) - B * atan((C*A)/(B*R));
                t2 = A/2 * log((R+B)/(R-B)) + B/2 * log((R+A)/(R-A)) - C * atan((A*B)/(C*R));
                coeff[0] += pow((-1),(double)(i+j+k)) * t0;
                coeff[1] += pow((-1),(double)(i+j+k)) * t1;
                coeff[2] += pow((-1),(double)(i+j+k)) * t2;
            }
        }
    }
    coeff[0] *= 1/(4*PI);
    coeff[1] *= 1/(4*PI);
    coeff[2] *= 1/(4*PI);

    return coeff;
}
