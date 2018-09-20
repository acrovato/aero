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

#define ANSI_COLOR_BLUE    "\x1b[1;34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;
using namespace Eigen;

int post(double sRef, double alpha, double Minf, Network &bPan, Field &fPan, double cL, double cD) {

    //// Begin post-processing
    cout << ANSI_COLOR_BLUE;
    cout << "******************************" << endl;
    cout << "*Beginning post-processing...*" << endl;
    cout << "******************************";
    cout << ANSI_COLOR_RESET << endl;

    // Paths
    string outPath;
    outPath = "../IO/";

    //// Write to file
    write_sp(outPath, sRef, alpha, Minf, bPan, cL, cD);
    write_fv(outPath, alpha, Minf, fPan);

    //// End post-processing
    cout << ANSI_COLOR_BLUE;
    cout << "*****************************" << endl;
    cout << "*Post-processing successful!*" << endl;
    cout << "*****************************";
    cout << ANSI_COLOR_RESET << endl;
    return 0;
}