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

//// Surface grid reading
// Read data from /IO/*.pts file which contain coordinates of points defining the corner of each panel
//
// I/O:
// - path: path to *.pts file
// - sGrid: temporary dynamic array to store panel vertices
// - bPan: body panel (structure)

#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include "read_sgrid.h"

#define NDIM 3

using namespace std;
using namespace Eigen;

void read_sgrid(string path, MatrixX3d &sGrid, Network &bPan){

    ifstream infile(path);
    string line;

    int i = 2;
    int a = 0;
    string section;
    bool secFLAG = 1;

    if(!infile.is_open())
    {
        cout << "File not found: " << path << endl;
        exit(EXIT_FAILURE);
    }
    else
        cout << "Found grid file: " << path << endl;

    getline(infile, line);
    cout << "Reading file... '" << line << "'" << endl;

    while (getline(infile, line))
    {

        if (secFLAG) {
            stringstream ss(line);
            ss >> section;
            secFLAG = 0;
        }

        else {

            if (section == "$size") {
                stringstream ss(line);
                ss >> bPan.nC >> bPan.nS;
                sGrid.resize(bPan.nC * bPan.nS, NDIM);
                secFLAG = 1;
            }

            else if (section == "$points") {
                stringstream ss(line);
                ss >> sGrid(a,0) >> sGrid(a,1) >> sGrid(a,2);
                a++;
            }

            else {
                cout << "Invalid section name: " << section << " at line " << i-1 << endl;
                exit(EXIT_FAILURE);
            }

        }
        i++;
    }

    // Set number of panels
    bPan.nC_ = bPan.nC - 1;
    bPan.nS_ = bPan.nS - 1;
    bPan.nP = bPan.nC_ * bPan.nS_;

    cout << "Done reading surface sGrid file!" << endl;
    cout << "Number of chordwise points: " << bPan.nC << endl;
    cout << "Number of spanwise points: " << bPan.nS << endl;
    cout << "Number of panels: " << bPan.nP << endl;
    #ifdef VERBOSE
        cout << "Surface points: " << sGrid.rows() << 'X' << sGrid.cols() << endl;
        for (int i = 0; i < bPan.nC*bPan.nS; ++i)
            cout << i << ' ' << sGrid(i,0) << ' ' << sGrid(i,1) << ' ' << sGrid(i,2) << endl;
    #endif
    cout << endl;
}