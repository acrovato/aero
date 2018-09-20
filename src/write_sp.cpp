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

//// Surface pressure writing
// Write panel pressure surface in external .dat and .pos files
//
// I/O:
// - outPath: path to write file
// - sRef: reference surface of the full wing
// - Minf: freestream Mach number
// - alpha: freestream angle of attack
// - bPan: (network of) body panels
// - cL: lift coefficient
// - cD: drag coefficient

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "write_sp.h"

#define  PI 3.14159

using namespace std;
using namespace Eigen;

void write_sp(string outPath, double sRef, double alpha, double Minf, Network &bPan, double cL, double cD) {

    //// Write to .dat file
    cout << "Writing surface pressure file in 'sp.dat'... " << flush;
    string outPathDat = outPath + "sp.dat";
    int idx = 0;
    ofstream sp;
    sp.open(outPathDat);

    // General information (header)
    sp << "Surface pressure file" << endl << endl;
    sp.width(25); sp << left << "Angle of attack: "; sp.width(20); sp << left << alpha*180/PI << endl;
    sp.width(25); sp << left << "Mach: "; sp.width(20); sp << left << Minf << endl;
    sp.width(25); sp << left << "Panels: "; sp.width(20); sp << left << bPan.nP << endl;
    sp.width(25); sp << left << "Full reference surface: "; sp.width(20); sp << left << sRef << endl;
    sp.width(25); sp << left << "Lift coefficient: "; sp.width(20); sp << left << cL << endl;
    sp.width(25); sp << left << "Drag coefficient: "; sp.width(20); sp << left << cD << endl;
    sp.width(25); sp << left << "Lift to drag ratio: "; sp.width(20); sp << left << cL/cD << endl << endl;

    // Panel pressure coefficient
    for (int j = 0; j < bPan.nS_; ++j) {
        // Labeling
        sp << "*** Stat" << j << " - y = "<< bPan.CG(idx + 1,1) << " ***" << endl;
        sp.width(15); sp << right << "x";
        sp.width(15); sp << right << "z";
        sp.width(15); sp << right << "cp" << endl;
        for (int i = 0; i < bPan.nC_; ++i) {
            idx = i + j * bPan.nC_;
            sp.width(15); sp << right << bPan.CG(idx,0);
            sp.width(15); sp << right << bPan.CG(idx,2);
            sp.width(15); sp << right << bPan.cP(idx) << endl;
        }
        sp << endl;
    }

    // Close file
    sp.close();
    cout << "Done!" << endl;

    //// Write to gmsh .pos file
    cout << "Writing surface pressure file in 'Cp.pos'... " << flush;
    string outPathPos = outPath + "Cp.pos";
    idx = 0;
    sp.open(outPathPos);

    // Gmsh header
    sp << "View \"Cp\" {" << endl;

    // Panel pressure coefficient
    for (int j = 0; j < bPan.nS_; ++j) {
        for (int i = 0; i < bPan.nC_; ++i) {
            sp << "SQ(";
            for (int k = 0; k < 3; ++k)
                sp << bPan.v0(idx,k) << ",";
            for (int k = 0; k < 3; ++k)
                sp << bPan.v1(idx,k) << ",";
            for (int k = 0; k < 3; ++k)
                sp << bPan.v2(idx,k) << ",";
            for (int k = 0; k < 2; ++k)
                sp << bPan.v3(idx,k) << ",";
            sp << bPan.v3(idx,2) << "){";
            for (int k = 0; k < 3; ++k)
                sp << bPan.cP(idx) << ",";
            sp << bPan.cP(idx) << "};" << endl;
            idx++;
        }
    }

    // Gmsh footer
    sp << "};" << endl << endl;

    // Close file
    sp.close();
    cout << "Done!" << endl;
}